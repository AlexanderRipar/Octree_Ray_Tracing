#include "och_string_util.h"

#include <string>

namespace och
{
	std::string pad(std::string s, std::string::size_type padded_size, bool align_right, bool cut_excess, char padding)
	{
		if (s.size() > padded_size)
			if (cut_excess)
				if (align_right)
					return s.substr(s.size() - padded_size);
				else
					return s.substr(0, padded_size);
			else
				return s;

		std::string total_padding;

		int sz = padded_size - s.size();

		for (int i = 0; i < sz; ++i)
			total_padding += padding;

		if (align_right)
			return total_padding + s;
		else
			return s + total_padding;
	}

	std::string abbreviate_byte_size(size_t bytes)
	{
		char size_magnitude;
		float factor;

		if (bytes < 1024)
		{
			size_magnitude = ' ';
			factor = 1.0F;
		}
		else if (bytes < 1024 * 1024)
		{
			size_magnitude = 'k';
			factor = 1.0F / 1024.0F;
		}
		else if (bytes < 1024 * 1024 * 1024)
		{
			size_magnitude = 'm';
			factor = 1.0F / (1024.0F * 1024.0F);
		}
		else
		{
			size_magnitude = 'g';
			factor = 1.0F / (1024.0F * 1024.0F * 1024.0F);
		}

		return std::to_string(static_cast<float>(bytes * factor)).substr(0, 6) + ' ' + size_magnitude + 'b';
	}

	int unchecked_hexval(char c)
	{
		return (c & 15) + ((c & 64) >> 6) * 9;
	}
}
