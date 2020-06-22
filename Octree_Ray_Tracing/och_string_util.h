#pragma once

#include <string>
#include <iostream>

static constexpr char nl = '\n';

namespace och
{
	std::string pad(std::string s, std::string::size_type padded_size, bool align_right = true, bool cut_excess = false, char padding = ' ');

	std::string abbreviate_byte_size(size_t bytes);

	//Assumes that c is a valid hex-character [0-9] [a-f] [A-F]
	int unchecked_hexval(char c);

	template<typename T>
	std::string int_to_binary_str(const T t)
	{
		std::string ret(sizeof(T) * 8, 'X');

		for (size_t i = 0; i < sizeof(T) * 8; ++i)
			ret[sizeof(T) * 8 - i - 1] = ('0' + ((t & (static_cast<T>(1) << i)) != 0));

		return ret;
	}

	template<char Zero, char One, typename T>
	std::string int_to_binary_str(const T t)
	{
		std::string ret(sizeof(T) * 8, 'X');

		for (size_t i = 0; i < sizeof(T) * 8; ++i)
			ret[sizeof(T) * 8 - i - 1] = t & (static_cast<T>(1) << i) ? One : Zero;

		return ret;
	}

	template<typename T>
	void bout(const T t)
	{
		std::cout << int_to_binary_str(t) << std::endl;
	}
}
