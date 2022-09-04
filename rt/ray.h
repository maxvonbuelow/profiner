#pragma once

#include <cmath>

struct Ray {
	float org[3];
	float dir[3];
	inline void normalize()
	{
		float n = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
		for (int i = 0; i < 3; ++i) dir[i] /= n;
	}
};
