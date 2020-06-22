#include "och_voxel.h"

#include <cstdint>
#include <cstdio>
#include <string>

namespace och
{
	int first_non_space(FILE* file)
	{
		int c = fgetc(file);

		while (isspace(c))
			c = fgetc(file);

		return c;
	}

	void process_name(FILE* file, int curr_vx, voxel_data& data)
	{
		int c = first_non_space(file);

		int i = 0;

		for (i = 0; i != name::max_len && c != ':'; ++i, c = fgetc(file))
			data.names[curr_vx].str[i] = c;

		if (i == 1)
		{
			printf("\nVoxel-names must contain at least one character\n");
			exit(0);
		}
		else if (i == name::max_len)
		{
			printf("\nVoxel-names may not exceed %i characters\n", name::max_len);
			exit(0);
		}

		for (; i < 16; ++i)
			data.names[curr_vx].str[i] = '\0';
	}

	void process_colour(FILE* file, int curr_vx, voxel_data& data)
	{
		for (int dir = 0; dir != 6; ++dir)
		{
			int c = first_non_space(file);

			for (int byte = 0; byte != 3; ++byte)
			{
				int val = unchecked_hexval(c) << 4;

				c = fgetc(file);

				val |= unchecked_hexval(c);

				c = fgetc(file);

				data.colours[curr_vx * 6 + dir].arr[byte] = val;
			}
		}
	}

	voxel_data read_voxel_data(const std::string& filename)
	{
		FILE* file;
		fpos_t file_beg;

		fopen_s(&file, filename.c_str(), "r");
		fgetpos(file, &file_beg);

		//Count number of voxels by counting ':'-occurrences
		int c, voxel_num = 0;

		do {
			c = fgetc(file);
			if (c == ':')
				++voxel_num;
		} while (c != EOF);

		voxel_data data(voxel_num);

		fsetpos(file, &file_beg);

		for (int curr_vx = 0; curr_vx != voxel_num; ++curr_vx)
		{
			process_name(file, curr_vx, data);

			process_colour(file, curr_vx, data);
		}

		return data;
	}
}