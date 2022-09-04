#pragma once

#include "ray.h"
#include <limits>

inline bool tri_intersect(float &t, float &uu, float &vv, const Ray &ray, const float *v0, const float *v1, const float *v2)
{
	// from wikipedia moeller trumbore
	const double EPSILON = std::numeric_limits<double>::epsilon();

	double e1x = v1[0] - v0[0], e1y = v1[1] - v0[1], e1z = v1[2] - v0[2];
	double e2x = v2[0] - v0[0], e2y = v2[1] - v0[1], e2z = v2[2] - v0[2];

	double hx = ray.dir[1] * e2z - ray.dir[2] * e2y;
	double hy = ray.dir[2] * e2x - ray.dir[0] * e2z;
	double hz = ray.dir[0] * e2y - ray.dir[1] * e2x;

	double a = e1x * hx + e1y * hy + e1z * hz;
	if (a > -EPSILON && a < EPSILON) return false; // This ray is parallel to this triangle.

	double f = 1.0 / a;
	double sx = ray.org[0] - v0[0], sy = ray.org[1] - v0[1], sz = ray.org[2] - v0[2];
	double u = f * (sx * hx + sy * hy + sz * hz);
	if (u < 0.0 || u > 1.0) return false;

	double qx = sy * e1z - sz * e1y;
	double qy = sz * e1x - sx * e1z;
	double qz = sx * e1y - sy * e1x;

	double v = f * (ray.dir[0] * qx + ray.dir[1] * qy + ray.dir[2] * qz);
	if (v < 0.0 || u + v > 1.0) return false;
	uu = u;
	vv = v;

	// At this stage we can compute t to find out where the intersection point is on the line.
	t = f * (e2x * qx + e2y * qy + e2z * qz);
	if (t > EPSILON) return true; // ray intersection
	else return false; // This means that there is a line intersection but not a ray intersection.
}
