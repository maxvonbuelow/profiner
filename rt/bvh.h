#pragma once
#include <cstdint>
#include <iostream>
#include <fstream>
#include "aabb.h"
#include <unordered_set>
#include <queue>
#include <list>
#include <stack>


enum Heuristic { SAH, MEDIAN };
enum Volume { V_AABB, V_SPHERE };

enum SplitDescent { NODE_LEFT, NODE_RIGHT };
struct BVHBuilder
{
	uint32_t maxlvl;
	void setmaxlvl(uint32_t ml)
	{
		maxlvl = ml;
		std::cout << "SETTING ML: " << maxlvl << std::endl;
	}
	uint32_t num_nodes()
	{
		return subtrees.size();
	}
	std::vector<AABB> bounds;
	std::vector<uint32_t> subtrees;
	struct Parent {
		uint32_t idx, subtree;
		SplitDescent desc;
		Parent() : subtree(0)
		{}
	};
	std::vector<Parent> parents;
	BVHBuilder() : parents(128) {}
	std::vector<uint32_t> leaf_nodes;
	uint32_t depth = 0;
	template <typename T>
	void serialize(std::ostream &os, const std::vector<T> &data)
	{
		uint64_t size = data.size();
		os.write((const char*)&size, 8);
		os.write((const char*)data.data(), data.size() * sizeof(T));
	}
	template <typename T>
	void unserialize(std::istream &is, std::vector<T> &data)
	{
		uint64_t size;
		is.read((char*)&size, 8);
		data.resize(size);
		is.read((char*)data.data(), data.size() * sizeof(T));
	}
	void backup(const char *name, int tag)
	{
		std::ofstream os((std::string("bvhbuilder") + name + "." + std::to_string(tag) + ".bin").c_str());
		std::vector<uint32_t> dummy;
		dummy.push_back(maxlvl);
		serialize<uint32_t>(os, dummy);
		serialize<AABB>(os, bounds);
		serialize<uint32_t>(os, subtrees);
		serialize<uint32_t>(os, leaf_nodes);
	}
	bool restore(const char *name, int tag)
	{
		std::cout << "RESTORE " << (std::string("bvhbuilder") + name + "." + std::to_string(tag) + ".bin").c_str() << std::endl;
		std::ifstream is((std::string("bvhbuilder") + name + "." + std::to_string(tag) + ".bin").c_str());
		if (!is) return false;
		std::vector<uint32_t> dummy;
		unserialize<uint32_t>(is, dummy);
		maxlvl = dummy.back();
		unserialize<AABB>(is, bounds);
		unserialize<uint32_t>(is, subtrees);
		unserialize<uint32_t>(is, leaf_nodes);
		return true;
	}
	uint32_t emit_node(uint32_t level = 0, uint32_t parent = -1u, SplitDescent desc = NODE_LEFT)
	{
		depth = std::max(level, depth);
		subtrees.emplace_back(0);
		uint32_t s = parents.size();
		while (level >= s) s <<= 1;
		parents.resize(s);
		parents[level].idx = parent;
		parents[level].subtree = 0;
		parents[level].desc = desc;

		while (level > 0) {
			++parents[level].subtree;
			if (parents[level].desc == NODE_LEFT) {
				subtrees[parents[level].idx] &= 0xc0000000u;
				subtrees[parents[level].idx] |= parents[level].subtree;
			}
			--level;
		}
		++parents[level].subtree;

		return num_nodes() - 1;
	}
	void set_axis(uint32_t nidx, uint32_t axis)
	{
		subtrees[nidx] |= axis << 30;
	}

	void set_bounds(uint32_t nidx, const AABB &l, const AABB &r)
	{
		bounds.push_back(l);
		bounds.push_back(r);
	}
	void set_leaf(uint32_t nidx, uint32_t *sat, uint32_t n, uint32_t nmax)
	{
		subtrees[nidx] = n;
		set_axis(nidx, 3); // mark as leaf


		uint32_t off = leaf_nodes.size();
		leaf_nodes.resize(off + nmax, -1u);
		std::copy(sat, sat + n, leaf_nodes.begin() + off);
	}

	void construct(float *cens, float *aabbs, float *verts, uint32_t n, uint32_t nleaf, Heuristic heuristic, Volume volume);

	void test(Volume vol, float *verts, uint32_t nmax)
	{
		std::cout << "Testing BVH" << std::endl;
		struct SE {
			uint32_t ni, li;
			std::list<uint32_t> bounds;
		};
		std::stack<SE> stack;

		uint32_t ni = 0, li = 0;
		stack.push(SE{ ni, li });
		do {
			SE t = stack.top(); stack.pop();
			ni = t.ni; li = t.li;

			uint32_t st = subtrees[ni];
			uint32_t axis = st >> 30, left_subtree = st & 0x3fffffffu;

			if (axis == 3) {
// 				std::cout << "sf" << std::endl;
				uint32_t nn = left_subtree;
				for (uint32_t b : t.bounds) {
					AABB aabb = bounds[b];
					for (int i = 0; i < nn; ++i) {
						uint32_t tri = leaf_nodes[li * nmax + i];
						for (int j = 0; j < 3; ++j) {
// 					std::cout << "sf.1 " << li << " " << i << " " << tri << std::endl;
							if (vol == V_AABB) {
								if (!aabb.test(verts + tri * 9 + j * 3)) {
									std::cout << "BVH ERROR tri " << tri << " not in " << b << std::endl;
									for (int jj = 0; jj < 3; ++jj) {
										for (int kk = 0; kk < 3; ++kk) {
											std::cout << (verts + tri * 9 + jj * 3)[kk] << " ";
										}
										std::cout << std::endl;
									}

									std::cout << "AABB:" << std::endl;
									for (int kk = 0; kk < 3; ++kk) {
										std::cout << aabb.min[kk] << " ";
									}
									std::cout << std::endl;
									for (int kk = 0; kk < 3; ++kk) {
										std::cout << aabb.max[kk] << " ";
									}
									std::cout << std::endl;
									std::exit(1);
								}
// 					std::cout << "sf.2" << std::endl;
							} else {
								if (!aabb.test_as_sphere(verts + tri * 9 + j * 3)) {
									std::cout << "BVH ERROR: (" << (verts + tri * 9 + j * 3)[0] << " " << (verts + tri * 9 + j * 3)[1] << " " << (verts + tri * 9 + j * 3)[1] << ") is not in (" << aabb.min[0] << " " << aabb.min[1] << " " << aabb.min[2] << " r=" << aabb.max[0] << ")" << std::endl;
									std::exit(1);
								}
							}
						}
					}
				}
// 				std::cout << "sf2" << std::endl;
			} else {
				uint32_t bi = ni - li;

				uint32_t cl = ni + 1, cr = ni + 1 + left_subtree;
				uint32_t ll = li, lr = li + (left_subtree + 1) / 2;

				stack.push(SE{ cr, lr, t.bounds });
				stack.top().bounds.push_back(bi * 2 + 1);
				stack.push(SE{ cl, ll, t.bounds });
				stack.top().bounds.push_back(bi * 2);
			}
		} while (!stack.empty());
		std::cout << "Testing BVH succ" << std::endl;
	}
};
