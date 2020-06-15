#pragma once

#include <cstdint>

namespace och
{
	constexpr uint64_t x_zero_mask_3d_64 = 0xEDB6DB6DB6DB6DB6LL;	//Sign bit set
	constexpr uint64_t y_zero_mask_3d_64 = 0xDB6DB6DB6DB6DB6DLL;
	constexpr uint64_t z_zero_mask_3d_64 = 0xB6DB6DB6DB6DB6DBLL;


	constexpr uint32_t x_zero_mask_3d_32 = 0xB6DB6DB6;
	constexpr uint32_t y_zero_mask_3d_32 = 0xEDB6DB6D;				//Sign bit set
	constexpr uint32_t z_zero_mask_3d_32 = 0xDB6DB6DB;


	constexpr uint16_t x_zero_mask_3d_16 = 0xEDB6;					//Sign bit set
	constexpr uint16_t y_zero_mask_3d_16 = 0xDB6D;
	constexpr uint16_t z_zero_mask_3d_16 = 0xB6DB;




	__forceinline uint32_t z_encode_8(uint8_t x, uint8_t y, uint8_t z);

	__forceinline uint64_t z_encode_16(uint16_t x, uint16_t y, uint16_t z);

	__forceinline uint64_t z_encode_32(uint32_t x, uint32_t y, uint32_t z);

	__forceinline uint16_t z_encode_8(uint8_t x, uint8_t y);

	__forceinline uint32_t z_encode_16(uint16_t x, uint16_t y);

	__forceinline uint64_t z_encode_32(uint32_t x, uint32_t y);
	




	uint32_t z_encode_8_noninline(uint8_t x, uint8_t y, uint8_t z);

	uint64_t z_encode_16_noninline(uint16_t x, uint16_t y, uint16_t z);

	uint64_t z_encode_32_noninline(uint32_t x, uint32_t y, uint32_t z);

	uint16_t z_encode_8_noninline(uint8_t x, uint8_t y);

	uint32_t z_encode_16_noninline(uint16_t x, uint16_t y);

	uint64_t z_encode_32_noninline(uint32_t x, uint32_t y);





	uint32_t z_splice_8_x(uint8_t x);

	uint32_t z_splice_8_y(uint8_t y);

	uint32_t z_splice_8_z(uint8_t z);


	uint64_t z_splice_16_x(uint16_t x);

	uint64_t z_splice_16_y(uint16_t y);

	uint64_t z_splice_16_z(uint16_t z);


	uint64_t z_splice_32_x(uint32_t x);

	uint64_t z_splice_32_y(uint32_t y);

	uint64_t z_splice_32_z(uint32_t z);
}