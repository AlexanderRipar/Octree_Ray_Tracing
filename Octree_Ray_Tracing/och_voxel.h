#pragma once

#include <cstdint>
#include <string>

#include "och_string_util.h"

//EXPECTED FILE-FORMAT
//
//voxel_name:
//colour for x_pos
//	-||-y_pos
//	-||-z_pos
//	-||-x_neg
//	-||-y_neg
//	-||-z_neg
//
//NOTES:
//-Any character between end of z_neg/file-start and ':' is treated as part of the Voxel-name, excluding leading whitespace.
//
//-Voxel-names may at most be 16 characters long (including whitespace and non-print)
//
//-Voxel-names must be unique and at least one character long
//
//- colours are to be specified as hex: RRGGBB. Alpha is automatically set to 255
//
//- indentation and whitespace is optional, as long as there is a whitespace character separating entries (not strictly necessary between z_neg and voxel_name)

namespace och
{
	struct pixel
	{
		union
		{
			uint8_t r, g, b, a;
			uint8_t arr[4]{ 0, 0, 0, 0xFF };
		};
	};

	struct name
	{
		static constexpr int max_len = 15;

		char str[max_len + 1];
	};

	struct voxel_data
	{
	private:

		const std::string filename;

		int voxel_cnt;
		pixel* colours;
		name* names;

	public:
		int get_cnt() const;

		const pixel* const get_colours() const;

		const name* get_names() const;

		void reload(std::string& errmsg);

		voxel_data(const std::string& filename, std::string& errmsg);

		~voxel_data();

	};
}
