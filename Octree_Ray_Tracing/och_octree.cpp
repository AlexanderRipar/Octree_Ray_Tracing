#include "och_octree.h"

#include "och_z_order.h"
#include "och_debug.h"

#include <immintrin.h>
#include <intrin.h>
#include <cstdio>



namespace och
{
	octree::octree(uint16_t depth, uint32_t table_capacity) : _table(create_table(table_capacity)), depth(depth), dim(1 << depth), table_capacity(table_capacity) { }

	octree::~octree()
	{
		delete[] _table;
	}

	octree::node* octree::create_table(int table_capacity)
	{
		node* table = new node[table_capacity];

		for (int i = 0; i != table_capacity; ++i)
			for (auto& c : table[i].children)
				c = 0;

		for (int i = 1; i != table_capacity; ++i) table[i].children[0] = i + 1;

		table[table_capacity - 1].children[0] = 0;

		return table;
	}

	bool octree::node::operator==(const node& n) const
	{
		return _mm256_movemask_epi8(_mm256_cmpeq_epi32(_mm256_load_si256(reinterpret_cast<const __m256i*>(children)), _mm256_load_si256(reinterpret_cast<const __m256i*>(n.children)))) == 0xFFFFFFFF;
	}

	bool octree::node::is_empty() const
	{
		return _mm256_movemask_epi8(_mm256_cmpeq_epi32(_mm256_load_si256(reinterpret_cast<const __m256i*>(children)), _mm256_setzero_si256())) == 0xFFFFFFFF;
	}

	uint32_t octree::alloc()
	{
		OCH_IF_DEBUG(++node_cnt);

		if (!head)
		{
			printf("\nToo many allocations...\n");
			exit(0);
		}

		uint32_t old_head = head;

		head = *reinterpret_cast<uint32_t*>(_table + head);

		for (auto& i : _table[old_head].children) i = 0;

		return old_head;
	}

	void octree::dealloc(uint32_t idx)
	{
		OCH_IF_DEBUG(--node_cnt);

		*reinterpret_cast<uint32_t*>(_table + idx) = head;

		head = idx;
	}

	void octree::set(int16_t x, int16_t y, int16_t z, uint32_t vx)
	{
		uint64_t index = och::z_encode_16(x, y, z);

		uint32_t curr = 0;//Root

		for (int i = depth - 1; i != 0; --i)
		{
			int node_idx = (index >> (3 * i)) & 7;

			if (!_table[curr].children[node_idx])
				_table[curr].children[node_idx] = alloc();

			curr = _table[curr].children[node_idx];
		}

		_table[curr].children[index & 7] = vx;
	}

	void octree::unset(int16_t x, int16_t y, int16_t z)
	{
		uint64_t index = och::z_encode_16(x, y, z);

		uint32_t curr = 0;//Root

		int sptr = 0;
		
		uint32_t stack[13];
		
		for (int i = depth - 1; i != 0; --i)
		{
			const int node_idx = (index >> (3 * i)) & 7;
		
			uint32_t child = _table[curr].children[node_idx];
		
			if (!child)
				return;
		
			stack[sptr] = curr;
		
			++sptr;

			curr = child;
		}
		
		_table[curr].children[index & 7] = 0;
		
		if (_table[curr].is_empty())
			dealloc(curr);
		else
			return;

		for (int i = 1; i != depth; ++i)
		{
			const int node_idx = (index >> (3 * i)) & 7;

			--sptr;

			_table[stack[sptr]].children[node_idx] = 0;

			if (_table[stack[sptr]].is_empty())
				dealloc(stack[sptr]);
			else
				return;
		}
	}

	uint32_t octree::at(int16_t x, int16_t y, int16_t z) const
	{
		uint64_t index = och::z_encode_16(x, y, z);

		uint32_t curr = 0;//Root

		for (int i = depth - 1; i > 0; --i)
		{
			int node_idx = (index >> (3 * i)) & 7;

			if (!_table[curr].children[node_idx])
			{
				return 0;
			}

			curr = _table[curr].children[node_idx];
		}

		return _table[curr].children[index & 7];
	}

	int octree::get_node_cnt() const
	{
		return node_cnt;
	}

	void octree::sse_trace(float ox, float oy, float oz, float dx, float dy, float dz, direction& hit_direction, uint32_t& hit_voxel, float& hit_time) const
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

		uint32_t node_idx = 0;															//Current node, initially set to root

		int level = 1;																	//Octree-"level"; Used as termination-criterion

		int min_t_idx = 8;																//Invalid initial index to help identify origin inside voxel

		__m128 _t_min = _mm_setzero_ps();

		//LABELLED LOOP

	PUSH:

		if (uint32_t child = _table[node_idx].children[(idx ^ inv_signs) & 7])
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

			hit_voxel = 0;

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

	void octree::sse_trace(och::float3 o, och::float3 d, direction& hit_direction, uint32_t& hit_voxel, float& hit_time) const
	{
		sse_trace(o.x, o.y, o.z, d.x, d.y, d.z, hit_direction, hit_voxel, hit_time);
	}

	//d := abs(rcp(d))
	//o := (d < 0.0F ? o : 3.0F - o) * d
	void octree::prepare_avx_trace(float* ox, float* oy, float* oz, float* dx, float* dy, float* dz, int cnt) const
	{
		
	}

	void octree::octree_prepare_avx_rest(float* ox, float* oy, float* oz, float* dx, float* dy, float* dz, int cnt) const
	{
		
	}

	void octree::avx_trace(float* ox, float* oy, float* oz, float* dx, float* dy, float* dz, direction* hit_direction, uint32_t* hit_voxel, float* hit_time, int cnt) const
	{

	}
}
