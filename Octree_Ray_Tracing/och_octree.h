#pragma once

#include <cstdint>

#include "och_tree_helper.h"
#include "och_vec.h"

namespace och
{
	class octree
	{
		//STRUCTS
	public:

		struct node
		{
			uint32_t children[8] alignas(32);

			bool operator==(const node& n) const;

			bool is_empty() const;
		};

		//DATA
	private:

		uint32_t head = 1;
		int node_cnt = 1;				//Take root into account

	public:

		node* const _table;
		const uint16_t depth;
		const uint16_t dim;
		const uint32_t table_capacity;

		//FUNCTIONS
	private:

		uint32_t alloc();

		void dealloc(uint32_t idx);

		node* create_table(int table_capacity);

	public:

		void set(int16_t x, int16_t y, int16_t z, uint32_t vx);

		void unset(int16_t x, int16_t y, int16_t z);

		uint32_t at(int16_t x, int16_t y, int16_t z) const;

		int get_node_cnt() const;
		
		void sse_trace(float ox, float oy, float oz, float dx, float dy, float dz, direction& hit_direction, uint32_t& hit_voxel, float& hit_time) const;

		void sse_trace(och::float3 o, och::float3 d, direction& hit_direction, uint32_t& hit_voxel, float& hit_time) const;

		void prepare_avx_trace(float* ox, float* oy, float* oz, float* dx, float* dy, float* dz, int cnt) const;

		void octree_prepare_avx_rest(float* ox, float* oy, float* oz, float* dx, float* dy, float* dz, int cnt) const;

		void avx_trace(float* ox, float* oy, float* oz, float* dx, float* dy, float* dz, direction* hit_direction, uint32_t* hit_voxel, float* hit_time, int cnt) const;

		octree(uint16_t depth, uint32_t table_capacity);

		~octree();
	};
}
