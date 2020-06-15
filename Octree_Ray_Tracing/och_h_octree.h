#pragma once

#include <cstdint>

#include "och_tree_helper.h"
#include "och_vec.h"

namespace och
{
	struct h_node_8
	{
		uint32_t children[8] alignas(32);

		bool operator==(const h_node_8& n) const;

		bool is_zero() const;
	};

	const h_node_8* get_node(uint32_t idx);

	void set_root(uint32_t idx);

	uint32_t register_h_node_8(h_node_8* node);

	void remove_h_node_8(uint32_t node_index);

	int get_node_count();

	int get_fill_count();

	int get_table_size();

	void clear_all();

	void sse_trace_h_octree(
		int depth, 
		float ox, 
		float oy, 
		float oz, 
		float dx, 
		float dy, 
		float dz, 
		direction& 
		hit_direction, 
		voxel& hit_voxel,
		float& hit_time);

	void sse_trace_h_octree(
		int depth, 
		och::float3 o, 
		och::float3 d, 
		direction& hit_direction, 
		voxel& hit_voxel,
		float& hit_time);

	void inv_trace_h_octree(
		int child_depth, 
		float ox, 
		float oy, 
		float oz, 
		float dx, 
		float dy, 
		float dz, 
		direction& hit_direction, 
		voxel& hit_voxel, 
		float& hit_time);

	void inv_trace_h_octree(
		int child_depth, 
		float3 o, 
		float3 d, 
		direction& hit_direction, 
		voxel& hit_voxel, 
		float& hit_time);
}
