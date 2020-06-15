#include "och_float.h"

#include <cstdint>

namespace och
{
	uint32_t float_as_u32(float f)
	{
		return *reinterpret_cast<uint32_t*>(&f);
	}

	int32_t float_as_i32(float f)
	{
		return *reinterpret_cast<int32_t*>(&f);
	}

	float u32_as_float(uint32_t f)
	{
		return *reinterpret_cast<float*>(&f);
	}

	float i32_as_float(int32_t f)
	{
		return *reinterpret_cast<float*>(&f);
	}

	std::string float_to_binary(float f, char zero, char one, char field_separator)
	{
		std::string s;

		uint32_t bitpat = float_as_u32(f);

		s.push_back(bitpat & (1 << 31) ? one : zero);

		s.push_back(field_separator);

		for (int i = 30; i >= 23; --i)
			s.push_back(bitpat & (1 << i) ? one : zero);

		s.push_back(field_separator);

		for (int i = 22; i >= 0; --i)
			s.push_back(bitpat & (1 << i) ? one : zero);

		return s;
	}

	std::string float_to_binary(uint32_t bitpat, char zero, char one, char field_separator)
	{
		std::string s;

		s.push_back(bitpat & (1 << 31) ? one : zero);

		s.push_back(field_separator);

		for (int i = 30; i >= 23; --i)
			s.push_back(bitpat & (1 << i) ? one : zero);

		s.push_back(field_separator);

		for (int i = 22; i >= 0; --i)
			s.push_back(bitpat & (1 << i) ? one : zero);

		return s;
	}

	float abs(float f)
	{
		uint32_t temp = float_as_u32(f) & (static_cast<uint32_t>(-1) >> 1);
		return *reinterpret_cast<float*>(&temp);
	}

	double abs(double f)
	{
		return *reinterpret_cast<uint64_t*>(&f) & (static_cast<uint64_t>(-1) >> 1);
	}

	long double abs(long double f)
	{
		return *reinterpret_cast<uint64_t*>(&f) & (static_cast<uint64_t>(-1) >> 1);
	}

	float sgn(float f)
	{
		*reinterpret_cast<uint32_t*>(&f) &= 1 << 31;
		*reinterpret_cast<uint32_t*>(&f) |= 0x3F80'0000;
		return f;
	}

	float sign_mask(float f)
	{
		return float_as_u32(f) & (1 << 31);
	}

	int8_t sgn_i8(float f)
	{
		return static_cast<int8_t>(float_as_i32(f) >> 31) | 1;
	}

	int16_t sgn_i16(float f)
	{
		return static_cast<int16_t>(float_as_i32(f) >> 31) | 1;
	}

	int32_t sgn_i32(float f)
	{
		return (float_as_i32(f) >> 31) | 1;
	}

	int64_t sgn_i64(float f)
	{
		return ((static_cast<int64_t>(*reinterpret_cast<int32_t*>(&f)) << 32) >> 63) | 1;
	}

	//Return true for negative numbers, and false for positive numbers
	bool sgn_bit(float f)
	{
		return float_as_u32(f) & (1 << 31);
	}

	//Returns the sign of f in the least significant bit (- : 1   //   + : 0)
	int32_t sgn_bit_i32(float f)
	{
		return static_cast<int32_t>((float_as_u32(f) & (1 << 31)) >> 31);
	}

	//Returns the sign of f in the most significant bit (- : 0x80000000   //   + : 0)
	int32_t sgn_bit_i32_pos31(float f)
	{
		return float_as_i32(f) & (1 << 31);
	}

	float min(float a, float b)
	{
		int x = float_as_i32(a);
		int y = float_as_i32(b);
		
		y = y + ((x - y) & ((x - y) >> 31));
		
		return *reinterpret_cast<float*>(&y);
	}

	float clear_mantissa(float f)
	{
		return i32_as_float(float_as_i32(f) & 0xFF800000);
	}

	bool is_in_upper_half(float f, int pivot)
	{
		uint32_t bitpat = float_as_u32(f + (1 << pivot));				//bitwise reinterpretation from input-float to uint32_t

		uint32_t mantissa = bitpat & 0x7FFFFFu;							//the mantissa-bits of f

		uint32_t relevant_exponent = (bitpat & 0x3F800000u) >> 23;		//the seven lower bits from the exponent of f

		return (mantissa >> (21 + pivot - relevant_exponent)) & 0x1;

		//return float_as_u32(f + pivot) & 0x800000U;
	}
}