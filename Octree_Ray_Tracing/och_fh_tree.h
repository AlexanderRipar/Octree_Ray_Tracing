#pragma once

#include <cstdint>
#include <immintrin.h>

#include "och_z_order.h"

namespace och
{
	constexpr int Log2_table_capacity = 19;
	constexpr int Depth = 12;
	constexpr int Log2_base_dim = 5;

	class fh_octree
	{
		//TEMPLATED CONSTANTS

	public:

		static constexpr int depth = Depth;
		static constexpr int dim = 1 << Depth;
		static constexpr int log2_table_capacity = Log2_table_capacity;
		static constexpr int table_capacity = 1 << Log2_table_capacity;
		static constexpr int base_dim = 1 << Log2_base_dim;

		static_assert(dim >= base_dim);
		static_assert(dim <= UINT16_MAX);

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

		struct hashtable
		{
			hashtable()
			{
				for (auto& i : cashes) i = 0;
				for (auto& i : refcounts) i = 0;
			}

			uint8_t cashes[table_capacity];

			uint32_t refcounts[table_capacity];

			node nodes[table_capacity];
		};

		struct basetable
		{
			uint32_t roots[base_dim * base_dim * base_dim];

			uint32_t& at(uint16_t x, uint16_t y, uint16_t z)
			{
				return roots[och::z_encode_16(x, y, z)];
			}
		};

	public:

		static constexpr size_t table_bytes = sizeof(hashtable);

		//DATA

		hashtable* const table = new hashtable;

		basetable* const base = new basetable;

		uint32_t fillcnt = 0;
		uint32_t nodecnt = 0;
		uint32_t max_refcnt = 0;
		
		//FUNCTIONS

		
	};
}