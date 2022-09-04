// Copyright (c) 2021, Max von Buelow, GRIS, Technical University of Darmstadt

#include "bvh.h"

static void sigtrap() { __asm__("int3"); }
#define __BP sigtrap()

#include <random>
#include <numeric>
#include <vector>
#include <algorithm>
#include <queue>
#include <stack>
#include <unordered_set>
#include <iostream>
#include "vec.h"
#include "bounding_sphere.h"

#define NBINS 256


struct SplitX {
	uint32_t o, l;
	uint32_t node;
	uint32_t split_axis;
	bool leaf;
	uint32_t level;
	SplitDescent desc;
	SplitX()
	{}
	SplitX(uint32_t o, uint32_t l, uint32_t node, uint32_t split_axis, bool leaf, uint32_t level = 0, SplitDescent desc = NODE_LEFT) : o(o), l(l), node(node), split_axis(split_axis), leaf(leaf), level(level), desc(desc)
	{}
};

void BVHBuilder::construct(float *cens, float *aabbs, float *verts, uint32_t n, uint32_t nleaf, Heuristic heuristic, Volume volume)
{
// 	std::cout << "CONSTRUCT BVH" << std::endl;
	int max_axis = 3; // TODO
	uint32_t leafminsplitcount = nleaf + (nleaf & 1) + 2;

// 	std::cout << "Overall surface: " << surface_all << std::endl;

	std::vector<uint32_t> perm(n);
	std::vector<uint32_t> tree;
	std::iota(perm.begin(), perm.begin() + n, 0);

// 	std::cout << "Initial compontents: " << compontents << std::endl;

	std::deque<SplitX> S;
	
	bool dosplit = n > nleaf;

	S.emplace_back(0, n, -1u, 0, !dosplit, 0);
	int o, l, nidx, split_axis, level, parent_node;

	double total = 0.;
	uint32_t finalized = 0;
	int maxlevel = 0;
	int RR = 0;
	int doPrintCnt = 0;
	while (!S.empty()) {
		if (doPrintCnt % 100 == 0) {
			std::cout << "\rStack size: " << S.size();
			std::cout.flush();
		}
		++doPrintCnt;
		SplitX s = S.back(); S.pop_back();
		o = s.o; l = s.l; parent_node = s.node; split_axis = s.split_axis; level = s.level;
		maxlevel = std::max(maxlevel, level);
		uint32_t cur_node = this->emit_node(level, parent_node, s.desc);
// 		std::cout << "len: " << l << std::endl;
// 		std::cout << "\rCur len: " << l << "                         ";
		std::cout.flush();

// 		// get bounding triangles and swap them to the first two positions of perm
		float val = std::numeric_limits<float>::infinity();
		uint32_t idx;
		// min
// 		std::cout << "SPLIT AXIS: " << split_axis << std::endl;
		if (s.leaf) {
			this->set_leaf(cur_node, perm.data() + o, l, nleaf);
			continue;
		}

// 		BIN1 | BIN2 | BIN3 | BIN4
// 		4 bins -> 3 splits
// 		
// 		l = 1 * binsizef
// 		l = 2 * binsizef
// 		l = 3 * binsizef

		static const int NSPLITS = NBINS - 1;
		uint32_t binsize = (l + NBINS - 1) / NBINS;
		float binsizef = (float)l / NBINS;
// 		if (binsizef < 1) binsizef = 1;
		int splits[NSPLITS];
		for (int j = 0; j < NSPLITS; ++j) {
			splits[j] = std::max((int)((j + 1) * binsizef), 1);
		}
		float min_cost = std::numeric_limits<float>::infinity();
		std::vector<uint32_t> dimperms(l * 3);
		float costs[NSPLITS * 3];
		vec3 mins_l[NSPLITS * 3], maxs_l[NSPLITS * 3], mins_r[NSPLITS * 3], maxs_r[NSPLITS * 3];
		float max_d = 0, max_d_axis;
		for (int axis = 0; axis < max_axis; ++axis) {
			for (int i = 0; i < l; ++i) {
				dimperms[i + axis * l] = perm[o + i];
			}

			uint32_t *dp = dimperms.data() + axis * l;
			std::sort(dp, dp + l, [&](uint32_t a, uint32_t b) {
				return cens[a * 3 + axis] < cens[b * 3 + axis];
			});

			float surface_l[NSPLITS] = { 0 }, surface_r[NSPLITS] = { 0 };
			float surface;

			if (volume == V_SPHERE) { // TODO: median
// 				std::cout << "SF1" << std::endl;
// 				int lastsplit = 0;
#pragma omp parallel for
				for (int j = 0; j < NSPLITS; ++j) {
					BoundingSphere sp_l(l);
// 					std::cout << "sf1" << std::endl;
					int split = splits[j];
					if (j && split == splits[j - 1]) {
// 					if ((j + 1) * binsize > l) {
// 						// no right-hand side
						surface_l[j] = std::numeric_limits<float>::infinity();
						surface_r[j] = std::numeric_limits<float>::infinity();
						continue;
					}
// 					lastsplit = split;
					for (int k = 0; k < split; ++k) {
						if (k >= l) {
							std::cout << "NBINS: " << NBINS << " l: " << l << std::endl;
							std::cerr << "Exit " << k << " " << l << " bin " << j << " binsize=" << binsize  << std::endl;
							std::exit(1); }
						sp_l.add3(verts + dp[k] * 9);
					}
// 					std::cout << "sf2" << std::endl;
					BoundingSphere sp_r(l);
					for (int k = split; k < l; ++k) {
						sp_r.add3(verts + dp[k] * 9);
					}
// 					if ((j + 1) * binsize >= l) { std::cout << (j + 1) * binsize << " " << l <<"|" << j << "/" << NBINS << std::endl; std::exit(1);}
// 					std::cout << "sf3" << std::endl;
					float sp_l_d[4], sp_r_d[4];
					sp_l.compute(sp_l_d + 1, sp_l_d);
// 					std::cout << "sf4" << std::endl;
					sp_r.compute(sp_r_d + 1, sp_r_d);
					surface_l[j] = sp_l_d[0] * sp_l_d[0];
					surface_r[j] = sp_r_d[0] * sp_r_d[0];
// 					if ((j + 1) * binsize >= l) surface_r[j] = std::numeric_limits<float>::infinity();
// 						if (std::isnan(surface_l[j / binsize]) || std::isinf(surface_l[j / binsize])) std::exit(1);
// 						if (std::isnan(surface_r[j / binsize]) || std::isinf(surface_r[j / binsize])) std::exit(1);

// 					std::cout << "sf5" << std::endl;
					mins_l[axis * NSPLITS + j] = vec3(sp_l_d[1], sp_l_d[2], sp_l_d[3]);
					maxs_l[axis * NSPLITS + j] = vec3(sp_l_d[0], 0.f, 0.f);
					mins_r[axis * NSPLITS + j] = vec3(sp_r_d[1], sp_r_d[2], sp_r_d[3]);
					maxs_r[axis * NSPLITS + j] = vec3(sp_r_d[0], 0.f, 0.f);
				}
// 				std::cout << "SF2" << std::endl;
				BoundingSphere sall(l);
				for (int j = 0; j < l; ++j) {
					sall.add3(verts + dp[j] * 9);
				}
// 				std::cout << "SF3" << std::endl;
				float spa_d[4];
				sall.compute(spa_d + 1, spa_d);
				surface = spa_d[0] * spa_d[0];
// 				std::cout << "SF4" << std::endl;
				if (spa_d[0] > max_d) {
					max_d = spa_d[0];
					max_d_axis = axis;
				}
			} else {
				vec3 min_l(std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity()), max_l(-std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity());
				vec3 min_r(std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity()), max_r(-std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity());

				int lastleft = 0;
				for (int j = 0; j < NSPLITS; ++j) {
					for (int k = lastleft; k < splits[j]; ++k) {
						min_l = min3(min_l, vec3(aabbs + dp[k] * 6));
						max_l = max3(max_l, vec3(aabbs + dp[k] * 6 + 3));
					}
					lastleft = splits[j];

					vec3 d = max_l - min_l;
					surface_l[j] = d.x * d.y + d.x * d.z + d.y * d.z;

					mins_l[axis * NSPLITS + j] = min_l;
					maxs_l[axis * NSPLITS + j] = max_l;
				}
				int lastright = l;
				for (int j = NSPLITS - 1; j >= 0; --j) {
					for (int k = splits[j]; k < lastright; ++k) {
						min_r = min3(min_r, vec3(aabbs + dp[k] * 6));
						max_r = max3(max_r, vec3(aabbs + dp[k] * 6 + 3));
					}
					lastright = splits[j];

					vec3 d = max_r - min_r;
					surface_r[j] = d.x * d.y + d.x * d.z + d.y * d.z;

					mins_r[axis * NSPLITS + j] = min_r;
					maxs_r[axis * NSPLITS + j] = max_r;
				}

// 				int left_bin = 0;
// 				int right_bin = NSPLITS - 1;
// 				for (int j = 0; j < binsize * NSPLITS; ++j) {
// 					int jj = binsize * NBINS - 1 - j;
// 
// 					if (j < l) {
// 						min_l = min3(min_l, vec3(aabbs + dp[j] * 6));
// 						max_l = max3(max_l, vec3(aabbs + dp[j] * 6 + 3));
// 					}
// 
// 					if (jj < l) {
// 						min_r = min3(min_r, vec3(aabbs + dp[jj] * 6));
// 						max_r = max3(max_r, vec3(aabbs + dp[jj] * 6 + 3));
// 					}
// 
// // 					if (jj == splits[right_bin]) {
// 					if (jj % binsize == 0) {
// // 						--right_bin;
// 						vec3 d = max_r - min_r;
// 						surface_r[jj / binsize - 1] = d.x * d.y + d.x * d.z + d.y * d.z;
// 
// 						mins_r[axis * NSPLITS + jj / binsize - 1] = min_r;
// 						maxs_r[axis * NSPLITS + jj / binsize - 1] = max_r;
// 					}
// // 					if (jj == splits[left_bin]) {
// // 						++left_bin;
// 					if (j % binsize == binsize - 1) {
// 						vec3 d = max_l - min_l;
// 						surface_l[j / binsize] = d.x * d.y + d.x * d.z + d.y * d.z;
// // 						if (std::isnan(surface_l[j / binsize]) || std::isinf(surface_l[j / binsize])) std::exit(1);/*surface_l[j / binsize] = std::numeric_limits<float>::infinity();*/
// 
// 						mins_l[axis * NSPLITS + j / binsize] = min_l;
// 						maxs_l[axis * NSPLITS + j / binsize] = max_l;
// 					}
// 				}
				vec3 d = max3(max_l, max_r) - min3(min_l, min_r);
				if (axis == 0) {
					if (d.x > max_d) {
						max_d = d.x;
						max_d_axis = 0;
					}
					if (d.y > max_d) {
						max_d = d.y;
						max_d_axis = 1;
					}
					if (d.z > max_d) {
						max_d = d.z;
						max_d_axis = 2;
					}
				}
				surface = d.x * d.y + d.x * d.z + d.y * d.z;
			}
			for (int k = 0; k < NSPLITS; ++k) {
				float lhs = surface_l[k] / surface;
				float rhs = surface_r[k] / surface;
// 				uint32_t num_l = (k + 1) * binsize;
				float num_l = splits[k];
				uint32_t num_r = l - num_l;
				costs[axis * NSPLITS + k] = lhs * num_l + rhs * num_r;
// 				costs[axis * NSPLITS + k] += (1 - std::min(num_l, leafminsplitcount) / (float)leafminsplitcount) + (1 - std::min(num_r, leafminsplitcount) / (float)leafminsplitcount);
			}
		}

		int best_idx = 0;
		for (int i = 0; i < NSPLITS * 3; ++i) {
			if (costs[i] < min_cost) {
				best_idx = i;
				min_cost = costs[i];
			}
		}
		int best_split = splits[(best_idx % NSPLITS)];
		int best_axis = best_idx / NSPLITS;
// 		std::cout << best_split << " " << best_axis << std::endl;

		if (heuristic == MEDIAN) {
			best_axis = max_d_axis/*(split_axis + 1) % 3*/;
			best_split = splits[NSPLITS / 2];
			best_idx = best_axis * NSPLITS + NSPLITS / 2;
		}



		for (int i = 0; i < l; ++i) {
			perm[o + i] = dimperms[i + best_axis * l];
		}

		// new split
		int nl = best_split;
		int nr = l - nl;

		// calculate complete AABBs
		AABB aabb_l(mins_l[best_idx].data(), maxs_l[best_idx].data()), aabb_r(mins_r[best_idx].data(), maxs_r[best_idx].data());
#if 0
		if (volume == V_SPHERE) {
// 			BoundingSphere sall(l);
// 			for (int j = 0; j < l; ++j) {
// 				sall.add3(verts + dp[j] * 9);
// 			}
// 			float spa_d[4];
// 			sall.compute(spa_d + 1, spa_d);
// 			surface = spa_d[0] * spa_d[0];
		} else {
			AABB X, Y;
			for (int i = 0; i < nl; ++i) {
				X.feed_min(aabbs + perm[o + i] * 6);
				X.feed_max(aabbs + perm[o + i] * 6 + 3);
			}
			for (int i = 0; i < nr; ++i) {
				Y.feed_min(aabbs + perm[o + nl + i] * 6);
				Y.feed_max(aabbs + perm[o + nl + i] * 6 + 3);
			}

			vec3 minl = mins_l[best_idx], maxl = maxs_l[best_idx];
			vec3 minr = mins_r[best_idx], maxr = maxs_r[best_idx];

// 	std::cout << "CONSTRUCT BVH" << std::endl;
			for (int i = 0; i < 3; ++i) {
				if (std::fabs(X.min[i] - minl[i]) > std::numeric_limits<float>::epsilon()) std::cerr << "UNEQ1 " << X.min[i] << " "<< minl[i]<< std::endl;
				if (std::fabs(X.max[i] - maxl[i]) > std::numeric_limits<float>::epsilon()) std::cerr << "UNEQ2" << std::endl;
				if (std::fabs(Y.min[i] - minr[i]) > std::numeric_limits<float>::epsilon()) std::cerr << "UNEQ3" << std::endl;
				if (std::fabs(Y.max[i] - maxr[i]) > std::numeric_limits<float>::epsilon()) std::cerr << "UNEQ4" << std::endl;
			}
		}
#endif
// 		if (volume == V_SPHERE) {
// 			// transform aabb into sphere
// 			vec3 min_l(aabb_l.min), min_r(aabb_r.min);
// 			vec3 dl = min_l - vec3(aabb_l.max);
// 			vec3 dr = min_r - vec3(aabb_r.max);
// 			vec3 cenl = dl / 2.f + min_l, cenr = dr / 2.f + min_r;
// 			float rl = std::max(dl.x, std::max(dl.y, dl.z));
// 			float rr = std::max(dr.x, std::max(dr.y, dr.z));
// 			for (int i = 0; i < 3; ++i) {
// 				aabb_l.min[i] = cenl[i];
// 				aabb_r.min[i] = cenr[i];
// 			}
// 			aabb_l.max[0] = rl;
// 			aabb_r.max[0] = rr;
// 		}


		uint32_t lhs_min, lhs_max, rhs_min, rhs_max;

		if (nl <= 0 || nr <= 0) {
			std::cout << std::endl << nl << " " << nr << std::endl;
			throw std::runtime_error("Negative split count");
		}
		this->set_axis(cur_node, best_axis);
		this->set_bounds(cur_node, aabb_l, aabb_r);

		bool r = nr > nleaf;
		bool l = nl > nleaf;

		S.emplace_back(o + nl, nr, cur_node, best_axis, !r, level + 1, NODE_RIGHT);
		S.emplace_back(o, nl, cur_node, best_axis, !l, level + 1, NODE_LEFT);
	}
	std::cout << std::endl;
// 	std::cout << "MAX LEVEL: " << maxlevel << std::endl;std::exit(1);
	this->setmaxlvl(maxlevel);
}
