#pragma once

#include <iostream>
#include <cstdint>

namespace och
{
	uint32_t float_as_u32(float f);

	int32_t float_as_i32(float f);

	float u32_as_float(uint32_t f);

	float i32_as_float(int32_t f);

	std::string float_to_binary(float f, char zero = '_', char one = 'X', char field_separator = '|');
	float abs(float f);

	std::string float_to_binary(uint32_t bitpat, char zero = '_', char one = 'X', char field_separator = '|');
	float abs(float f);

	double abs(double f);

	long double abs(long double f);

	float sgn(float f);

	int8_t sgn_i8(float f);

	int16_t sgn_i16(float f);

	int32_t sgn_i32(float f);

	int64_t sgn_i64(float f);

	bool sgn_bit(float f);

	int32_t sgn_bit_i32(float f);

	int32_t sgn_bit_i32_pos31(float f);

	float min(float a, float b);

	float clear_mantissa(float f);

	bool is_in_upper_half(float f, int pivot);

	//Returns true if a float within [0, 2 * Pivot[ is greater than Pivot; Otherwise, returns false
	template<uint32_t Pivot> bool is_in_upper_half_template(float f)
	{
		constexpr float bias = static_cast<float>(Pivot);

		f += bias;

		return *reinterpret_cast<uint32_t*>(&f) & 0x800000U;
	}

	struct accumulate_float_bits
	{
		void insert(float f)
		{
			uint32_t bits = float_as_u32(f);

			_union |= bits;
			_intersect &= bits;
		}

		uint32_t get_union() const
		{
			return _union;
		}

		uint32_t get_intersect() const
		{
			return _intersect;
		}

		uint32_t get_identifiers() const
		{
			return _union & _intersect;
		}

		uint32_t get_unique(accumulate_float_bits& o)
		{
			return (get_intersect() ^ o.get_intersect()) & (get_union() ^ o.get_union());
		}

		void reset()
		{
			_union = 0;
			_intersect = 0xFFFFFFFF;
		}

	private:

		uint32_t _union = 0;
		uint32_t _intersect = 0xFFFFFFFF;
	};

	template<uint32_t Pivot, uint32_t Granularity> int test_upper_half_on_range(float lo, float hi)
	{
		float range = hi - lo;

		bool res_union = false;
		bool res_intersect = true;

		for (int i = 0; i < Granularity; ++i)
		{
			float f = lo + ((float) i / (float) Granularity) * range;

			bool res = is_in_upper_half<Pivot>(f);

			res_union |= res;
			res_intersect &= res;
		}

		return res_union != res_intersect ? -1 : res_union;
	}

	//The idea here is that the bit exp bits into a float's mantissa is zero if it is in the lower half 
	//(i.e. in    [n * 2^ModPowTwo,     n * 2^ModPowTwo + 2 ^ (ModPowTwo - 1)[    )
	//and one otherwise
	template<uint32_t ModPowTwo> int is_in_upper_mod_half(float f)
	{
		static_assert(ModPowTwo != 0);

		constexpr float f_offset = static_cast<float>(1 << ModPowTwo);	//offset applied to f, as the algorithm does not work for f in [0, 2^ModPowTwo[ otherwise

		constexpr uint32_t pow_bias = 21 + ModPowTwo;					//precomputed shifting constant depending on "modulo" specified at compile-time

		uint32_t bitpat = float_as_u32(f + f_offset);					//bitwise reinterpretation from input-float to uint32_t

		uint32_t mantissa = bitpat & 0x7FFFFFu;							//the mantissa-bits of f

		uint32_t relevant_exponent = (bitpat & 0x3F800000u) >> 23;		//the seven lower bits from the exponent of f

		return (mantissa >> (pow_bias - relevant_exponent)) & 0x1;
	}
}