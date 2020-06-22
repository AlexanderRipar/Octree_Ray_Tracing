#pragma once

#include <cstdlib>
#include <cstdint>
#include <cstdio>
#include <immintrin.h>
#include <intrin.h>

#include "och_tree_helper.h"
#include "och_z_order.h"
#include "och_vec.h"
#include "och_debug.h"

namespace och
{

	template<int Log2_table_capacity, int Depth>
	class h_octree
	{
		//TEMPLATED CONSTANTS

	public:

		static constexpr int depth = Depth;
		static constexpr int dim = 1 << Depth;
		static constexpr int log2_table_capacity = Log2_table_capacity;
		static constexpr int table_capacity = 1 << Log2_table_capacity;
		static constexpr float voxel_dim = 1.0F / dim;

	private:

		static constexpr int idx_mask = ((table_capacity - 1) >> 4) << 4;

		//STRUCTS

	public:

		struct node
		{
			uint32_t children[8] alignas(32);

			bool operator==(const node& n) const
			{
				return _mm256_movemask_epi8(_mm256_cmpeq_epi32(_mm256_load_si256(reinterpret_cast<const __m256i*>(children)), _mm256_load_si256(reinterpret_cast<const __m256i*>(n.children)))) == 0xFFFFFFFF;
			}

			bool is_zero() const
			{
				return _mm256_movemask_epi8(_mm256_cmpeq_epi32(_mm256_load_si256(reinterpret_cast<const __m256i*>(children)), _mm256_setzero_si256())) == 0xFFFFFFFF;
			}

			uint32_t hash() const
			{
				constexpr uint32_t Prime = 0x01000193; //   16777619
				constexpr uint32_t Seed = 0x811C9DC5;  // 2166136261

				const char* node_bytes = reinterpret_cast<const char*>(children);

				uint32_t hash = Seed;

				for (int i = 0; i < 32; ++i)
					hash = (*node_bytes++ ^ hash) * Prime;

				return hash;
			}
		};

	private:

		struct node_hashtable
		{
			node_hashtable()
			{
				for (auto& i : cashes) i = 0;
				for (auto& i : refcounts) i = 0;
			}

			uint8_t cashes[table_capacity];

			uint32_t refcounts[table_capacity];

			node nodes[table_capacity];
		};

	public:

		static constexpr size_t table_bytes = sizeof(node_hashtable);

		//DATA

	private:

		node_hashtable* const table = new node_hashtable;

		uint32_t root_idx = 0;

		uint32_t fillcnt = 0;
		uint32_t nodecnt = 0;
		uint32_t max_refcnt = 0;

		//FUNCTIONS

	public:

		~h_octree()
		{
			delete[] table;
		}

		uint32_t register_node(const node& n)
		{
			if (fillcnt > static_cast<uint32_t>((table_capacity * 0.9375F)))
			{
				printf("\ntabled nodes: %i\n=> Table too full. Exiting...%c\n", fillcnt, 7);
				exit(0);
			}

			uint32_t hash = n.hash();

			uint32_t index = hash & idx_mask;

			uint8_t cash = static_cast<uint8_t>(hash >> log2_table_capacity);

			if (cash == 0)//Avoid empty
				cash = 1;
			else if (cash == 0xFF)//Avoid Gravestone
				cash = 0x7F;

			uint32_t last_grave = -1;

			while (table->cashes[index])
			{
				if (table->cashes[index] == 0xFF)
					last_grave = index;
				if (table->cashes[index] == cash)
					if (table->nodes[index] == n)
					{
						OCH_IF_DEBUG(++nodecnt);

						++table->refcounts[index];

						return index + 1;
					}

				index = (index + 1) & (table_capacity - 1);
			}

			OCH_IF_DEBUG(++nodecnt; ++fillcnt);

			if (last_grave != -1)
				index = last_grave;

			table->cashes[index] = cash;

			table->nodes[index] = n;

			table->refcounts[index] = 1;

			return index + 1;
		}

		void remove_node(const uint32_t idx)
		{
			--table->refcounts[idx - 1];

			OCH_IF_DEBUG(--nodecnt);

			if (!table->refcounts[idx - 1])
			{
				OCH_IF_DEBUG(--fillcnt);

				table->cashes[idx - 1] = 0xFF;//Gravestone
			}
		}

		void set(uint16_t x, uint16_t y, uint16_t z, uint32_t v)
		{
			if ((x | y | z) >= dim)
				return;

			uint64_t index = och::z_encode_16(x, y, z);

			uint32_t stk[depth];

			int d = depth - 1;

			//Build stack of nonzero nodes, starting from root, terminating at leaf/null
			for (uint32_t curr = root_idx; curr && d >= 0; --d)
			{
				int c_idx = (index >> (3 * d)) & 7;

				stk[d] = curr;

				curr = table->nodes[curr - 1].children[c_idx];
			}

			uint32_t child = v;

			int _d = 0;

			//If the stack did not reach leaf: starting at leaf, work up to (excluding) last nonnull node and insert them into table
			if (++d)
			{
				if (!v)
					return;

				while (_d != d)
				{
					node n{ 0, 0, 0, 0, 0, 0, 0, 0 };

					int c_idx = (index >> (3 * _d++)) & 7;

					n.children[c_idx] = child;

					child = register_node(n);
				}
			}

			//Delete, modify, and reinsert nodes from stack
			for (int i = d; i != depth; ++i)
			{
				remove_node(stk[i]);

				node n = table->nodes[stk[i] - 1];

				int c_idx = (index >> (3 * i)) & 7;

				n.children[c_idx] = child;

				if (n.is_zero())
					child = 0;
				else
					child = register_node(n);
			}

			set_root(child);
		}

		uint32_t at(int x, int y, int z)
		{
			uint64_t index = och::z_encode_16(x, y, z);

			uint32_t curr = root_idx;//Root

			for (int i = depth - 1; i != 0; --i)
			{
				int node_idx = (index >> (3 * i)) & 7;

				if (!table->nodes[curr - 1].children[node_idx])
				{
					return 0;
				}

				curr = table->nodes[curr - 1].children[node_idx];
			}

			return table->nodes[curr - 1].children[index & 7];
		}

		void set_root(uint32_t idx)
		{
			root_idx = idx;
		}

		uint32_t get_root()
		{
			return root_idx;
		}

		uint32_t get_fillcnt() const
		{
			return fillcnt;
		}

		uint32_t get_nodecnt() const
		{
			return nodecnt;
		}

		uint32_t get_max_refcnt() const
		{
			return max_refcnt;
		}

		void clear()
		{
			for (auto& i : table->cashes) i = 0;
		}

		//TRACING

		void sse_trace(float ox, float oy, float oz, float dx, float dy, float dz, direction& hit_direction, uint32_t& hit_voxel, float& hit_time) const
		{
			const __m128 _signbit = _mm_set1_ps(-0.0F);										//Numeric constants

			const __m128 _one_point_five = _mm_set1_ps(1.5F);

			const __m128 _three = _mm_add_ps(_one_point_five, _one_point_five);

			const __m128 _x_mask = _mm_castsi128_ps(_mm_set_epi32(0, 0, 0, -1));

			const __m128 _y_mask = _mm_castsi128_ps(_mm_set_epi32(0, 0, -1, 0));

			const __m128 _z_mask = _mm_castsi128_ps(_mm_set_epi32(0, -1, 0, 0));

			__m128 _d = _mm_set_ps(0.0F, dz, dy, dx);										//Populate simd-register with direction

			__m128 _o = _mm_set_ps(0.0F, oz, oy, ox);										//Populate simd-register with origin

			const __m128 _sign_mask = _mm_cmplt_ps(_mm_setzero_ps(), _d);					//Create sign mask from _d: Negative -> 0; Positive -> -1;

			_d = _mm_or_ps(_d, _signbit);													//Ensure that direction is negative

			_o = _mm_andnot_ps(_signbit, _mm_sub_ps(_mm_and_ps(_sign_mask, _three), _o));	//If direction was initially positive -> o = 3 - o (reflect around 1.5)

			const __m128 _coef = _mm_rcp_ps(_d);											//Get reciprocal of direction as ray-coefficient

			const __m128 _bias = _mm_xor_ps(_signbit, _mm_mul_ps(_coef, _o));				//Get -o * coef as ray-bias

			__m128 _pos = _mm_and_ps(_o, _one_point_five);									//Extract aligned position from origin

			const int inv_signs = _mm_movemask_ps(_sign_mask);								//Mask of direction-sign: d < 0 ? 0 : 1

			int idx = _mm_movemask_ps(_mm_cmpeq_ps(_pos, _one_point_five));					//Initial idx, derived from comparing pos and tree-middle

			__m128 _dim_bit = _mm_castsi128_ps(_mm_set1_epi32(1 << 22));					//Bit used for masking / adding / subtracting. It corresponds to the mantissa-bit for child-dimension

			uint32_t parents[depth - 1];													//Stack of parent-voxels for restoring on pop

			uint32_t* curr_paren = parents;													//Pointer to current parent in parent-stack

			uint32_t node_idx = root_idx;													//Current node, initially set to root

			int level = 1;																	//Octree-"level"; Used as termination-criterion

			int min_t_idx = 8;																//Invalid initial index to help identify origin inside voxel

			__m128 _t_min = _mm_setzero_ps();

			//LABELLED LOOP

		PUSH:

			if (uint32_t child = table->nodes[node_idx - 1].children[(idx ^ inv_signs) & 7])
			{
				if (level++ == depth)//HIT
				{
					hit_voxel = child;

					hit_direction = static_cast<direction>((min_t_idx >> 1) + 3 * ((inv_signs & min_t_idx) == 0));

					hit_time = _t_min.m128_f32[0];

					return;
				}

				*(curr_paren++) = node_idx;

				node_idx = child;

				_dim_bit = _mm_castsi128_ps(_mm_srai_epi32(_mm_castps_si128(_dim_bit), 1));

				__m128 _mid_node_pos = _mm_or_ps(_pos, _dim_bit);

				__m128 _t_at_mid_node = _mm_fmadd_ps(_mid_node_pos, _coef, _bias);

				__m128 _new_idx_mask = _mm_cmpge_ps(_t_at_mid_node, _t_min);

				idx = _mm_movemask_ps(_new_idx_mask);

				__m128 _pos_incr = _mm_and_ps(_dim_bit, _new_idx_mask);

				_pos = _mm_or_ps(_pos, _pos_incr);

				goto PUSH;
			}

		STEP:

			const __m128 _t = _mm_fmadd_ps(_pos, _coef, _bias);								//Calculate times at position

			__m128 _min_mask;

			unsigned int tx = _t.m128_u32[0];
			unsigned int ty = _t.m128_u32[1];
			unsigned int tz = _t.m128_u32[2];

			if (tx <= ty && tx <= tz)
			{
				min_t_idx = 1;
				_t_min = _mm_castsi128_ps(_mm_set1_epi32(tx));
				_min_mask = _x_mask;

			}
			else if (ty < tx && ty <= tz)
			{
				min_t_idx = 2;
				_t_min = _mm_castsi128_ps(_mm_set1_epi32(ty));
				_min_mask = _y_mask;
			}
			else
			{
				min_t_idx = 4;
				_t_min = _mm_castsi128_ps(_mm_set1_epi32(tz));
				_min_mask = _z_mask;
			}

			__m128 _pos_decr;																//Declared ahead of goto to shut up intellisense; Should not (significantly) influence generated code

			if (!(idx & min_t_idx))
				goto POP;

			_pos_decr = _mm_and_ps(_min_mask, _dim_bit);

			_pos = _mm_andnot_ps(_pos_decr, _pos);

			idx ^= min_t_idx;

			goto PUSH;

		POP:

			if (--level == 0)//MISS
			{
				hit_direction = direction::exit;

				hit_voxel = 0;

				hit_time = INFINITY;

				return;
			}

			node_idx = *(--curr_paren);														//Restore parent-node from parent-stack

			_pos = _mm_andnot_ps(_dim_bit, _pos);											//Restore parent-position

			_dim_bit = _mm_castsi128_ps(_mm_slli_epi32(_mm_castps_si128(_dim_bit), 1));		//Update dimension-bit

			const __m128 _pos_lowbit = _mm_and_ps(_dim_bit, _pos);

			const __m128 _pos_isupper = _mm_cmpeq_ps(_dim_bit, _pos_lowbit);

			idx = _mm_movemask_ps(_pos_isupper);											//Restore index

			goto STEP;
		}

		void sse_trace(och::float3 o, och::float3 d, direction& hit_direction, uint32_t& hit_voxel, float& hit_time) const
		{
			sse_trace(o.x, o.y, o.z, d.x, d.y, d.z, hit_direction, hit_voxel, hit_time);
		}

		/*void avx_trace(const float* ox, const float* oy, const float* oz, const float* dx, const float* dy, const float* dz, direction* hit_direction, uint32_t* hit_voxel, float* hit_time) const
		{
			const __m256 _zero = _mm256_setzero_ps();
			const __m256 _sbit = _mm256_set1_ps(-0.0F);
			const __m256 _three = _mm256_set1_ps(3.0F);
			const __m256 _one_point_five = _mm256_set1_ps(1.5F);

			__m256 _ox = _mm256_load_ps(ox);
			__m256 _oy = _mm256_load_ps(oy);
			__m256 _oz = _mm256_load_ps(oz);

			__m256 _dx = _mm256_load_ps(dx);//direction
			__m256 _dy = _mm256_load_ps(dy);
			__m256 _dz = _mm256_load_ps(dz);

			__m256 _sx = _mm256_cmp_ps(dx, _zero, 29);//sign
			__m256 _sy = _mm256_cmp_ps(dz, _zero, 29);
			__m256 _sz = _mm256_cmp_ps(dy, _zero, 29);

			_dx = _mm256_or_ps(_dx, _sbit);
			_dy = _mm256_or_ps(_dy, _sbit);
			_dz = _mm256_or_ps(_dz, _sbit);

			_ox = _mm_or_ps(_sbit, _mm_sub_ps(_mm_and_ps(_sx, _three), _ox));//origin
			_oy = _mm_or_ps(_sbit, _mm_sub_ps(_mm_and_ps(_sy, _three), _oy));
			_oz = _mm_or_ps(_sbit, _mm_sub_ps(_mm_and_ps(_sz, _three), _oz));

			const __m256 _cx = _mm256_rcp_ps(_dx);//coeff
			const __m256 _cy = _mm256_rcp_ps(_dy);
			const __m256 _cz = _mm256_rcp_ps(_dz);

			const __m256 _bx = _mm256_mul_ps(_cx, _ox);//bias
			const __m256 _by = _mm256_mul_ps(_cy, _oy);
			const __m256 _bz = _mm256_mul_ps(_cz, _oz);

			__m256 _px = _mm256_and_ps(_ox, _one_point_five);//position
			__m256 _py = _mm256_and_ps(_oy, _one_point_five);
			__m256 _pz = _mm256_and_ps(_oz, _one_point_five);

			int idx_x = _mm256_movemask_ps(_mm256_cmp_ps(_px, _one_point_five, 0);
			int idx_y = _mm256_movemask_ps(_mm256_cmp_ps(_py, _one_point_five, 0);
			int idx_z = _mm256_movemask_ps(_mm256_cmp_ps(_pz, _one_point_five, 0);

			int idx = och::z_encode_8(signs_x, signs_y, signs_z);
		}*/
	};
}
