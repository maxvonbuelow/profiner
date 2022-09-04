#pragma once

// A CGAL wrapper that may be replaced with something else to reduce dependencies...

//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Min_sphere_of_points_d_traits_3.h>
//#include <CGAL/Min_sphere_of_spheres_d.h> 

//typedef CGAL::Simple_cartesian<double>           K;
//typedef CGAL::Min_sphere_of_points_d_traits_3<K,double> Traits;
//typedef CGAL::Min_sphere_of_spheres_d<Traits>    Min_sphere;
//typedef K::Point_3                               Point;

struct BoundingSphere {
//	std::vector<Point> points;
	BoundingSphere(uint32_t pts)
	{
//		points.reserve(pts);
	}

	void add3(const float *threepts)
	{
		add(threepts);
		add(threepts + 3);
		add(threepts + 6);
	}
	void add(const float *pt)
	{
//		points.emplace_back(pt[0], pt[1], pt[2]);
	}
	void compute(float *cen, float *r)
	{
#if 0
		Min_sphere ms(points.data(), points.data() + points.size());
		Min_sphere::Cartesian_const_iterator ccib = ms.center_cartesian_begin();
		for (int i = 0; i < 3; ++i) {
			cen[i] = *ccib++;
		}
		*r = ms.radius() + /*10*std::numeric_limits<float>::epsilon()*/0.00001;
		for (auto pt : points) {
			float dx = pt[0] - cen[0], dy = pt[1] - cen[1], dz = pt[2] - cen[2];
			float dd = dx * dx + dy * dy + dz * dz;
			if (dd > *r * *r) {
				std::cout << "Point " << pt[0] << " " << pt[1] << " " << pt[2] << " outside " << cen[0] << " " << cen[1] << " " << cen[2] << " r=" << *r << " dist = " << std::sqrt(dd) << std::endl;
				std::exit(1);
			}
		}
#endif
	}
};
