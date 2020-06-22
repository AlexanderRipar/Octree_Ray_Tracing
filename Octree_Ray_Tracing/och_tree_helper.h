#pragma once

#include <cstdint>

namespace och
{
	enum class direction
	{
		x_pos = 0,
		y_pos = 1,
		z_pos = 2,
		x_neg = 3,
		y_neg = 4,
		z_neg = 5,
		exit = 6,
		inside = 7,
		error = 8
	};
}
