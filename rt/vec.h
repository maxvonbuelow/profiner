#pragma once

#include <ostream>

struct vec3 {
	float x, y, z;
	inline vec3() {}
	inline vec3(float a, float b, float c) : x(a), y(b), z(c) {}
	inline vec3(const float *a) : x(a[0]), y(a[1]), z(a[2]) {}
	inline vec3 operator-(const vec3 &other) { return vec3(x - other.x, y - other.y, z - other.z); }
	inline vec3 operator+(const vec3 &other) { return vec3(x + other.x, y + other.y, z + other.z); }
	inline vec3 operator/(float a) { return vec3(x / a, y / a, z / a); }
	inline vec3 operator*(float a) { return vec3(x * a, y * a, z * a); }
	inline float operator[](int a) const { return (&x)[a]; }
	inline float &operator[](int a) { return (&x)[a]; }
	inline float *data() { return &x; }
};
inline vec3 cross(const vec3 &a, const vec3 &b)
{
	return vec3(
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x
	);
}
inline float dot(const vec3 &a, const vec3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline float sign(float x) { return x == 0 ? 0 : x < 0 ? -1 : 1; }
inline float dot2(const vec3 &a) { return dot(a, a); }
inline float clamp(float v, float a, float b) { return std::min(std::max(v, a), b); }
// float abs(float x) { return x < 0 ? -x : x; }

inline vec3 min3(const vec3 &a, const vec3 &b) { return vec3(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)); }
inline vec3 max3(const vec3 &a, const vec3 &b) { return vec3(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)); }

inline float max3v(float a, float b, float c) { return std::max(std::max(a, b), c); }
inline float min3v(float a, float b, float c) { return std::min(std::min(a, b), c); }

inline float max2v(float a, float b) { return std::max(a, b); }
inline float min2v(float a, float b) { return std::min(a, b); } 

inline vec3 abs3(const vec3 &a) { return vec3(std::abs(a.x), std::abs(a.y), std::abs(a.z)); }
inline float min3i(const vec3 &a) { return min3v(a.x, a.y, a.z); }
inline float max3i(const vec3 &a) { return max3v(a.x, a.y, a.z); }

inline std::ostream &operator<<(std::ostream &os, const vec3 &v)
{
	os << '[' << v.x << ", " << v.y << ", " << v.z << ']';
	return os;
}
