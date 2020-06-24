#include "och_voxel.h"

#include <cstdint>
#include <cstdio>
#include <string>

namespace och
{
	int first_non_space(FILE* file)
	{
		int c;

		do
		{
			c = fgetc(file);

			if (c == EOF)
				return c;
		}
		while (isspace(c));

		return c;
	}

	int get_voxel_cnt(const std::string& filename, std::string& errmsg)
	{
		FILE* file;
		if (int file_error = fopen_s(&file, filename.c_str(), "r"))
		{
			char msg[128];

			if (!strerror_s(msg, file_error))
				errmsg = "Could not open " + filename + ": " + std::string(msg) + " (errno: " + std::to_string(file_error) + ")";
			else
				errmsg = "holy pooper, wahapen?";

			return 0;
		}

		//Count number of voxels by counting ':'-occurrences
		int c, voxel_cnt = 0;

		do {
			c = fgetc(file);

			if (c == ':')
				++voxel_cnt;
		} while (c != EOF);

		if (voxel_cnt == 0)
		{
			errmsg = "File did not contain a valid voxel";
		}

		fclose(file);

		return voxel_cnt;
	}

	int voxel_data::get_cnt() const
	{
		return voxel_cnt;
	}

	const pixel* const voxel_data::get_colours() const
	{
		return colours;
	}

	const name* voxel_data::get_names() const
	{
		return names;
	}

	void voxel_data::reload(std::string& errmsg)
	{
		int new_voxel_cnt = get_voxel_cnt(filename, errmsg);

		if (errmsg != "")//Set from get_voxel_cnt
			return;

		if (voxel_cnt != new_voxel_cnt)
		{
			errmsg = "New Voxel-count does not match old";
			return;
		}

		FILE* file;

		if (int file_error = fopen_s(&file, filename.c_str(), "r"))
		{
			char msg[128];

			if (!strerror_s(msg, file_error))
				errmsg = "Could not open " + filename + ": " + std::string(msg) + " (errno: " + std::to_string(file_error) + ")";
			else
				errmsg = "holy pooper, wahapen? Your error was seemingly too big (>= 128 characters). (errno: " + std::to_string(file_error) + ")";

			return;
		}

		for (int curr_vx = 0; curr_vx != voxel_cnt; ++curr_vx)
		{
			//Read name
			int c = first_non_space(file);

			int i = 0;

			for (i = 0; i != name::max_len && c != ':'; ++i, c = fgetc(file))
			{
				if (c == EOF)
				{
					errmsg = filename + " ended unexpectedly";
					fclose(file);
					return;
				}

				if (iscntrl(c))
				{
					errmsg = "Voxel-names may not contain control-characters (see voxel number" + std::to_string(curr_vx + 1) + ")";
				}

				names[curr_vx].str[i] = c;
			}

			if (i == 1)
			{
				errmsg = "Voxel-names must contain at least one character";
				fclose(file);
				return;
			}
			else if (i == name::max_len)
			{
				errmsg = "Voxel-names may not exceed " + std::to_string(name::max_len) + " characters";
				fclose(file);
				return;
			}

			for (; i < 16; ++i)
				names[curr_vx].str[i] = '\0';

			//Read colours
			for (int dir = 0; dir != 6; ++dir)
			{
				int c = first_non_space(file);

				for (int byte = 0; byte != 3; ++byte)
				{
					if (c == EOF)
					{
						errmsg = filename + " ended unexpectedly";
						fclose(file);
						return;
					}

					if (!isxdigit(c))
					{
						errmsg = "Non-hex character in colour-value (" + std::string(names[curr_vx].str) + " at colour no. " + std::to_string(dir + 1) + ")";
						fclose(file);
						return;
					}

					int val = unchecked_hexval(c) << 4;

					c = fgetc(file);

					if (c == EOF)
					{
						errmsg = filename + " ended unexpectedly";
						fclose(file);
						return;
					}

					if (!isxdigit(c))
					{
						errmsg = "Non-hex character in colour-value (" + std::string(names[curr_vx].str) + " at colour no. " + std::to_string(dir + 1) + ")";
						fclose(file);
						return;
					}

					val |= unchecked_hexval(c);

					c = fgetc(file);

					colours[curr_vx * 6 + dir].arr[byte] = val;
				}
			}
		}

		fclose(file);

		return;
	}

	voxel_data::voxel_data(const std::string& filename, std::string& errmsg) : filename(filename), voxel_cnt( get_voxel_cnt(filename, errmsg) ), colours( new pixel[6 * voxel_cnt] ), names( new name[voxel_cnt] )
	{
		FILE* file;

		if (errmsg != "")//Set from get_voxel_cnt()
			return;

		if (int file_error = fopen_s(&file, filename.c_str(), "r"))
		{
			char msg[128];
			
			if (!strerror_s(msg, file_error))
				errmsg = "Could not open " + filename + ": " + std::string(msg) + " (errno: " + std::to_string(file_error) + ")";
			else
				errmsg = "holy pooper, wahapen? Your error was seemingly too big (>= 128 characters). (errno: " + std::to_string(file_error) + ")";

			return;
		}

		for (int curr_vx = 0; curr_vx != voxel_cnt; ++curr_vx)
		{
			//Read name
			int c = first_non_space(file);

			int i = 0;

			for (i = 0; i != name::max_len && c != ':'; ++i, c = fgetc(file))
			{
				if (c == EOF)
				{
					errmsg = filename + " ended unexpectedly";
					fclose(file);
					return;
				}

				if (iscntrl(c))
				{
					errmsg = "Voxel-names may not contain control-characters (see voxel number" + std::to_string(curr_vx + 1) + ")";
				}

				names[curr_vx].str[i] = c;
			}

			if (i == 1)
			{
				errmsg = "Voxel-names must contain at least one character";
				fclose(file);
				return;
			}
			else if (i == name::max_len)
			{
				errmsg = "Voxel-names may not exceed " + std::to_string(name::max_len) + " characters";
				fclose(file);
				return;
			}

			for (; i < 16; ++i)
				names[curr_vx].str[i] = '\0';

			//Read colours
			for (int dir = 0; dir != 6; ++dir)
			{
				int c = first_non_space(file);

				for (int byte = 0; byte != 3; ++byte)
				{
					if (c == EOF)
					{
						errmsg = filename + " ended unexpectedly";
						fclose(file);
						return;
					}

					if (!isxdigit(c))
					{
						errmsg = "Non-hex character in colour-value (" + std::string(names[curr_vx].str) + " at colour no. " + std::to_string(dir + 1) + ")";
						fclose(file);
						return;
					}

					int val = unchecked_hexval(c) << 4;

					c = fgetc(file);

					if (c == EOF)
					{
						errmsg = filename + " ended unexpectedly";
						fclose(file);
						return;
					}

					if (!isxdigit(c))
					{
						errmsg = "Non-hex character in colour-value (" + std::string(names[curr_vx].str) + " at colour no. " + std::to_string(dir + 1) + ")";
						fclose(file);
						return;
					}

					val |= unchecked_hexval(c);

					c = fgetc(file);

					colours[curr_vx * 6 + dir].arr[byte] = val;
				}
			}
		}

		fclose(file);

		return;
	}

	voxel_data::~voxel_data()
	{
		delete[] colours;
		delete[] names;
	}
}
