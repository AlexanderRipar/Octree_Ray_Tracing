#include "och_h_octree.h"

#include <cstdlib>
#include <cstdio>
#include <immintrin.h>
#include <intrin.h>
#include <limits>

#include "och_tree_helper.h"
#include "och_vec.h"
#include "och_debug.h"

namespace och
{
	static constexpr int log2_tablesize = 24;

	static constexpr int tablesize = 1 << log2_tablesize;

	static constexpr int idx_mask = ((tablesize - 1) >> 4) << 4;

	uint32_t fill_count = 0;
	uint32_t node_count = 0;
	uint16_t max_refcount = 0;

	bool h_node_8::operator==(const h_node_8& n) const
	{
		return _mm256_movemask_epi8(_mm256_cmpeq_epi32(_mm256_load_si256(reinterpret_cast<const __m256i*>(children)), _mm256_load_si256(reinterpret_cast<const __m256i*>(n.children)))) == 0xFFFFFFFF;
	}

	bool h_node_8::is_zero() const
	{
		return _mm256_movemask_epi8(_mm256_cmpeq_epi32(_mm256_load_si256(reinterpret_cast<const __m256i*>(children)), _mm256_setzero_si256())) == 0xFFFFFFFF;
	}

	uint32_t hash_node(const h_node_8* node)
	{
		const uint32_t Prime = 0x01000193; //   16777619
		const uint32_t Seed = 0x811C9DC5; // 2166136261

		const char* node_bytes = reinterpret_cast<const char*>(node);

		uint32_t hash = Seed;

		for (int i = 0; i < 32; ++i)
			hash = (*node_bytes++ ^ hash) * Prime;

		return hash;
	}

	//DO NOT CONSTRUCT ON STACK IF YOU DON'T DESIRE DEATH AND MISERY AND PAIN AND SUCH! This thing is pretty chunky
	struct node_hashtable
	{
		uint8_t cashes[tablesize];

		uint16_t refcounts[tablesize];//Not sufficient for depth > 9

		h_node_8 nodes[tablesize];
	};

	node_hashtable table;

	uint32_t root_idx = 0;

	const h_node_8* get_node(uint32_t idx)
	{
		return table.nodes - 1 + idx;
	}

	void set_root(uint32_t idx)
	{
		root_idx = idx;
	}

	uint32_t register_h_node_8(h_node_8* node)
	{
		if (fill_count > static_cast<int>((tablesize * 0.9375F)))
		{
			printf("\ntabled nodes: %i\n=> Table too full. Exiting...\n", fill_count);
			exit(0);
		}

		uint32_t hash = hash_node(node);

		uint32_t index = hash & idx_mask;

		uint8_t cash = static_cast<uint8_t>((hash >> log2_tablesize) | 1);

		while (table.cashes[index])
		{
			if (table.cashes[index] == cash)
				if (table.nodes[index] == *node)
				{
					OCH_IF_DEBUG(++node_count);

					++table.refcounts[index];

					return index + 1;
				}

			index = (index + 1) & (tablesize - 1);
		}
		OCH_IF_DEBUG(++node_count; ++fill_count;)

			table.cashes[index] = cash;

		table.nodes[index] = *node;

		table.refcounts[index] = 1;

		return index + 1;
	}

	void remove_h_node_8(uint32_t node_index)
	{
		int refcount = --table.refcounts[node_index - 1];

		OCH_IF_DEBUG(--node_count;)

			if (!refcount)
			{
				OCH_IF_DEBUG(--fill_count;)

					table.cashes[node_index - 1] = 0;
			}
	}

	void clear_all()
	{
		for (int i = 0; i < tablesize; ++i)
			table.cashes[i] = 0;
	}

	int get_node_count()
	{
		return node_count;
	}

	int get_fill_count()
	{
		return fill_count;
	}

	int get_table_size()
	{
		return tablesize;
	}




	//Helper functions for tracing
	__forceinline direction get_hit_direction_for_inv(int sign_mask, int last_step)
	{
		int step_sign = sign_mask & last_step;

		switch (last_step)
		{
		case 0:
			return direction::inside;
		case 1:
			return step_sign ? direction::x_pos : direction::x_neg;
		case 2:
			return step_sign ? direction::y_pos : direction::y_neg;
		case 4:
			return step_sign ? direction::z_pos : direction::z_neg;
		default:
			return direction::error;
		}
	}

	__forceinline int get_last_step_for_inv(float last_time, float tx_pos, float ty_pos, float tx_coeff, float ty_coeff, float tx_bias, float ty_bias)
	{
		if (last_time == 0.0F)
			return 0;

		float tx = (tx_pos + 1.0F) * tx_coeff + tx_bias;
		float ty = (ty_pos + 1.0F) * ty_coeff + ty_bias;

		if (tx == last_time)
			return 1;
		else if (ty == last_time)
			return 2;
		else
			return 4;
	}







	//Tracing functions
	void avx_trace_h_tree(int depth, const float* ox, const float* oy, const float* oz, const float* dx, const float* dy, const float* dz, direction* out_direction, voxel* out_voxel)
	{
		const __m256 _signbit = _mm256_set1_ps(-0.0F);									//Numeric constants

		const __m256 _one_point_five = _mm256_set1_ps(1.5F);

		const __m256 _three = _mm256_add_ps(_one_point_five, _one_point_five);

		__m256 _dx = _mm256_load_ps(dx);
		__m256 _dy = _mm256_load_ps(dy);
		__m256 _dz = _mm256_load_ps(dz);

		__m256 _ox = _mm256_load_ps(ox);
		__m256 _oy = _mm256_load_ps(oy);
		__m256 _oz = _mm256_load_ps(oz);

		const __m256 _sx = _mm256_cmp_ps(_mm256_setzero_ps(), _dx, 17);					//Inverse signmasks; 0 < d ? 0xFFFF'FFFF : 0x0000'0000;
		const __m256 _sy = _mm256_cmp_ps(_mm256_setzero_ps(), _dy, 17);					
		const __m256 _sz = _mm256_cmp_ps(_mm256_setzero_ps(), _dz, 17);					

		_dx = _mm256_or_ps(_dx, _signbit);												
		_dy = _mm256_or_ps(_dy, _signbit);
		_dz = _mm256_or_ps(_dz, _signbit);

		_ox = _mm256_andnot_ps(_signbit, _mm256_sub_ps(_mm256_and_ps(_sx, _three), _ox));
		_oy = _mm256_andnot_ps(_signbit, _mm256_sub_ps(_mm256_and_ps(_sy, _three), _oy));
		_oz = _mm256_andnot_ps(_signbit, _mm256_sub_ps(_mm256_and_ps(_sz, _three), _oz));

		const __m256 _cx = _mm256_rcp_ps(_dx);											//Step-coefficients
		const __m256 _cy = _mm256_rcp_ps(_dy);
		const __m256 _cz = _mm256_rcp_ps(_dz);

		const __m256 _bx = _mm256_andnot_ps(_signbit, _mm256_mul_ps(_ox, _cx));			//Step-biases
		const __m256 _by = _mm256_andnot_ps(_signbit, _mm256_mul_ps(_oy, _cy));
		const __m256 _bz = _mm256_andnot_ps(_signbit, _mm256_mul_ps(_oz, _cz));

		__m256 _px = _mm256_and_ps(_ox, _one_point_five);								//Ray-positions
		__m256 _py = _mm256_and_ps(_oy, _one_point_five);
		__m256 _pz = _mm256_and_ps(_oz, _one_point_five);

		__m256 _dim_bits = _mm256_castsi256_ps(_mm256_set1_epi32(1 << 22));

		uint32_t parents[8][13];

		uint8_t curr_paren[8]{ 0, 0, 0, 0, 0, 0, 0, 0 };
	}



	/*
	void pck_trace_h_tree(uint32_t root_idx, int depth, const float* ox, const float* oy, const float* oz, const float* dx, const float* dy, const float* dz, int elem_cnt, direction* out_direction, voxel* out_voxel)
	{

	}
	*/


	void sse_trace_h_octree(int depth, float ox, float oy, float oz, float dx, float dy, float dz, direction& hit_direction, voxel& hit_voxel, float& hit_time)
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

		uint32_t parents[13];															//Stack of parent-voxels for restoring on pop

		uint32_t* curr_paren = parents;													//Pointer to current parent in parent-stack

		uint32_t node_idx = root_idx;														//Current node, initially set to root

		int level = 1;																	//Octree-"level"; Used as termination-criterion

		int min_t_idx = 8;																//Invalid initial index to help identify origin inside voxel

		__m128 _t_min = _mm_setzero_ps();

		//LABELLED LOOP

	PUSH:

		if (uint32_t child = get_node(node_idx)->children[(idx ^ inv_signs) & 7])
		{
			if (level++ == depth)//HIT
			{
				hit_voxel = static_cast<voxel>(child);

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

		if (!(idx & min_t_idx))
			goto POP;

		const __m128 _pos_decr = _mm_and_ps(_min_mask, _dim_bit);

		_pos = _mm_andnot_ps(_pos_decr, _pos);

		idx ^= min_t_idx;
		
		goto PUSH;

	POP:

		if (--level == 0)//MISS
		{
			hit_direction = direction::exit;

			hit_voxel = voxel::air;

			hit_time = 0.0F;

			return;
		}

		node_idx = *(--curr_paren);															//Restore parent-node from parent-stack

		_pos = _mm_andnot_ps(_dim_bit, _pos);											//Restore parent-position

		_dim_bit = _mm_castsi128_ps(_mm_slli_epi32(_mm_castps_si128(_dim_bit), 1));		//Update dimension-bit

		const __m128 _pos_lowbit = _mm_and_ps(_dim_bit, _pos);

		const __m128 _pos_isupper = _mm_cmpeq_ps(_dim_bit, _pos_lowbit);

		idx = _mm_movemask_ps(_pos_isupper);											//Restore index

		goto STEP;
	}

	void sse_trace_h_octree(int depth, och::float3 o, och::float3 d, direction& hit_direction, voxel& hit_voxel, float& hit_time)
	{
		sse_trace_h_octree(depth, o.x, o.y, o.z, d.x, d.y, d.z, hit_direction, hit_voxel, hit_time);
	}


	
	void inv_trace_h_octree(int child_depth, float ox, float oy, float oz, float dx, float dy, float dz, direction& hit_direction, voxel& hit_voxel, float& hit_time)
	{
		constexpr float inf = std::numeric_limits<float>::infinity();

		uint32_t parents[11];

		uint32_t* curr_paren = parents;

		uint32_t node_idx = root_idx;

		const int signbit_x = dx < 0.0F ? 0 : 1;//!och::sgn_bit(dx);//Inverted
		const int signbit_y = dy < 0.0F ? 0 : 1;//!och::sgn_bit(dy);
		const int signbit_z = dz < 0.0F ? 0 : 1;//!och::sgn_bit(dz);

		const int sign_mask = signbit_x | (signbit_y << 1) | (signbit_z << 2);//Also inverted

		float child_dim = static_cast<float>(1 << child_depth);

		const float max_dim = child_dim;

		if (signbit_x) dx = -dx, ox = 2 * max_dim - ox;
		if (signbit_y) dy = -dy, oy = 2 * max_dim - oy;
		if (signbit_z) dz = -dz, oz = 2 * max_dim - oz;

		const int idx_x = ox >= child_dim;
		const int idx_y = oy >= child_dim;
		const int idx_z = oz >= child_dim;

		int64_t idx = idx_x | (idx_y << 1) | (idx_z << 2);//Since this is a 32-bit int, the maximum child-dim is 9 ((9 + 1) * 3 == 30 Bits)

		const float tx_coeff = dx == 0 ? 0 : 1.0F / dx;
		const float ty_coeff = dy == 0 ? 0 : 1.0F / dy;
		const float tz_coeff = dz == 0 ? 0 : 1.0F / dz;

		const float tx_bias = tx_coeff == 0 ? inf : -ox * tx_coeff;
		const float ty_bias = ty_coeff == 0 ? inf : -oy * ty_coeff;
		const float tz_bias = tz_coeff == 0 ? inf : -oz * tz_coeff;

		float tx_pos = child_dim * idx_x;
		float ty_pos = child_dim * idx_y;
		float tz_pos = child_dim * idx_z;

		float last_time = 0.0F;//Minimum tw_pos

		//LABELLED LOOP

	PUSH:

		if (uint32_t child = get_node(node_idx)->children[(idx ^ sign_mask) & 7])//PUSH
		{
			if (child_dim == 1)//HIT
			{
				//if (!last_step) goto STEP;

				hit_voxel = static_cast<voxel>(child);

				hit_direction = get_hit_direction_for_inv(sign_mask, get_last_step_for_inv(last_time, tx_pos, ty_pos, tx_coeff, ty_coeff, tx_bias, ty_bias));

				hit_time = last_time;

				return;
			}

			*(curr_paren++) = node_idx;

			node_idx = child;

			child_dim *= 0.5F;

			bool n_idx_x = dx * last_time + ox >= tx_pos + child_dim;
			bool n_idx_y = dy * last_time + oy >= ty_pos + child_dim;
			bool n_idx_z = dz * last_time + oz >= tz_pos + child_dim;

			idx <<= 3;

			idx |= n_idx_x | (n_idx_y << 1) | (n_idx_z << 2);

			if (n_idx_x) tx_pos += child_dim;
			if (n_idx_y) ty_pos += child_dim;
			if (n_idx_z) tz_pos += child_dim;

			goto PUSH;
		}

	STEP:

		float tx = tx_pos * tx_coeff + tx_bias;
		float ty = ty_pos * ty_coeff + ty_bias;
		float tz = tz_pos * tz_coeff + tz_bias;

		if (tx <= ty && tx <= tz)
		{
			if (!(idx & 1))
				goto POP;
			else
				idx ^= 1;

			last_time = tx;
			tx_pos -= child_dim;
		}
		else if (ty < tx && ty <= tz)
		{
			if (!(idx & 2))
				goto POP;
			else
				idx ^= 2;

			last_time = ty;
			ty_pos -= child_dim;
		}
		else
		{
			if (!(idx & 4))
				goto POP;
			else
				idx ^= 4;

			last_time = tz;
			tz_pos -= child_dim;
		}

		goto PUSH;

	POP:

		if (child_dim == max_dim)//MISS
		{
			hit_direction = direction::exit;

			hit_voxel = voxel::air;

			hit_time = 0;

			return;
		}

		node_idx = *(--curr_paren);

		//TODO: Take care of tw_pos //Should be doable branchless with bit-ops
		if (idx & 1) tx_pos -= child_dim;
		if (idx & 2) ty_pos -= child_dim;
		if (idx & 4) tz_pos -= child_dim;

		child_dim *= 2;

		idx >>= 3;

		goto STEP;
	}

	void inv_trace_h_octree(int child_depth, float3 o, float3 d, direction& hit_direction, voxel& hit_voxel, float& hit_time)
	{
		inv_trace_h_octree(child_depth, o.x, o.y, o.z, d.x, d.y, d.z, hit_direction, hit_voxel, hit_time);
	}
}
