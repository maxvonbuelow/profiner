/*
 * Copyright (C) 2022, Max von Buelow
 * Technical University of Darmstadt - Interactive Graphics Systems Group
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#include "meminf_storage.h"
#include <cmath>
#include <atomic>
#include <iostream>
#include <vector>
#include <map>
#include "caches.h"
#include "cachesim_alloc.h"

#define ADDRSIZE 1
#define L1LINE 32
#define L2LINE 64
#define PERCENTIFY(p) ((int)(p * 10000) / 100.f)

enum MATag { L1HIT = 1, L2HIT = 2, ATOMIC = 4, IS_WRITE = 8, WRITE_TYPE_THROUGH = 16 };

struct MA { // 64 bit now!!!!
	uint32_t off;
	uint16_t bufid;
	uint16_t tags;

	void set_atomic(bool a)
	{
		tags &= ~(decltype(tags))ATOMIC;
		tags |= a ? ATOMIC : 0;
	}
	void set_l1hit(bool h)
	{
		tags &= ~(decltype(tags))L1HIT;
		tags |= h ? L1HIT : 0;
	}
	void set_l2hit(bool h)
	{
		tags &= ~(decltype(tags))L2HIT;
		tags |= h ? L2HIT : 0;
	}
	void set_write_through(bool write_through)
	{
		tags |= IS_WRITE;
		tags &= ~(decltype(tags))WRITE_TYPE_THROUGH;
		tags |= write_through ? WRITE_TYPE_THROUGH : 0;
	}
	bool goes_to_l2() const
	{
		return !get_l1hit() || is_writethrough();
	}
	void set_sm(int smid)
	{
		tags &= 31;
		tags |= smid << 5;
	}
	int get_sm() const
	{
		return tags >> 5;
	}

	bool get_atomic() const
	{
		return tags & ATOMIC;
	}
	bool get_l1hit() const
	{
		return tags & L1HIT;
	}
	bool get_l2hit() const
	{
		return tags & L2HIT;
	}
	bool is_write() const
	{
		return tags & IS_WRITE;
	}
	bool is_writethrough() const
	{
		return is_write() && (tags & WRITE_TYPE_THROUGH);
	}
};

struct MACoalescingInfo {
	uint32_t lanes;
	std::bitset<160> suboffs; // offsets in the 32 bytes (log2(32) * 32)

	uint8_t get_suboff(int i) const
	{
		uint8_t out = 0;
		int mask = 1;
		for (int j = 0; j < 5; ++j) {
			out |= suboffs[i + j] ? mask : 0;
			mask <<= 1;
		}
		return out;
	}
	void set_suboff(int i, uint8_t val)
	{
		int mask = 1;
		for (int j = 0; j < 5; ++j) {
			suboffs[i + j] = !!(val & mask);
			mask <<= 1;
		}
	}
};

struct Alloc {
	uint64_t size;
	int order;

	int sizehash_minor;
	uint64_t ld_reqs = 0;
	uint64_t st_reqs = 0;
	uint64_t ld_txn = 0;
	uint64_t st_txn = 0;
	std::vector<uint32_t> memaccs;

	uint32_t sh()
	{
		return (uint32_t)size << 8 | sizehash_minor;
	}
};



struct CacheConf {
	uint32_t l1_capacity;
	uint32_t l2_capacity;
	uint32_t l2_assoc;
};

struct CUDADim {
	uint32_t bx, by, bz;
	uint32_t tx, ty, tz;
};

struct CacheSim {
// 	typedef RandomEvictionCacheSim CacheType;
	typedef PLRUCacheSim CacheType;
	typedef NaiveRDCacheSim CacheL2Type;
// 	typedef LRUL2CacheSim CacheL2Type;
	std::vector<MA> memaccs;
	std::vector<uint64_t> gwarpids;
	std::vector<MACoalescingInfo> coalescing;
	std::vector<CacheType> reuse_l1_storage;
	CUDADim griddim;

	std::size_t size()
	{
		return memaccs.size();
	}
	int dims_at(uint64_t idx, CUDADim *dims)
	{
		CUDADim dim;
		uint64_t g_warp_id = gwarpids[idx];
		uint32_t lanes = coalescing[idx].lanes;

		// decompose global warp ID into block ID and local warp ID
		uint32_t n_threads = griddim.tx * griddim.ty * griddim.tz;
		uint32_t n_warps = (n_threads + 31) / 32;
		uint32_t block_id = g_warp_id / n_warps;
		uint32_t local_warp_id = g_warp_id % n_warps;

		// decompose block ID
		dim.bz = block_id / (griddim.bx * griddim.by);
		uint32_t tmp = block_id % (griddim.bx * griddim.by);
		dim.by = tmp / griddim.by;
		dim.bz = tmp % griddim.by;

		// calculate local thread IDs from local warp ID
		int j = 0;
		for (int i = 0; i < 32; ++i) {
			if (!(lanes & ((uint32_t)1 << i))) continue; // this memory reference is not included in that lane

			uint32_t l_thread_id = local_warp_id * 32 + i;

			// decompose local thread ID
			dim.tz = l_thread_id / (griddim.tx * griddim.ty);
			tmp = l_thread_id % (griddim.tx * griddim.ty);
			dim.ty = tmp / griddim.ty;
			dim.tx = tmp % griddim.ty;

			dims[j++] = dim;
		}
		return j;
	}
	void grid(uint32_t nbx, uint32_t nby, uint32_t nbz, uint32_t ntx, uint32_t nty, uint32_t ntz)
	{
		griddim = CUDADim{ nbx, nby, nbz,    ntx, nty, ntz };
	}

	uint32_t nwarps()
	{
		uint32_t nb = griddim.bx * griddim.by * griddim.bz;
		uint32_t nt = griddim.tx * griddim.ty * griddim.tz;
		uint32_t nw = (nt + 31) / 32;
		return nb * nw;
	}
	
	struct Conf {
		bool record_coalescing = false, record_gwarpid = false;
	};
	Conf conf;

	CacheConf cacheconf;

	CacheL2Type reuse_l2;
	void clear()
	{
		memaccs.clear();
	}
	void operator()(uint64_t off, uint8_t footprint, uint32_t bufid, uint8_t sm, const uint32_t *lanes, const uint8_t *offs, uint32_t nlanes/*, uint64_t allocoff*/, bool atomic, bool is_write, bool is_writethrough, uint64_t gwarpid)
	{
		memaccs.push_back(MA{ (uint32_t)(off / ADDRSIZE), (uint16_t)bufid });
		memaccs.back().set_sm(sm);
		memaccs.back().set_atomic(atomic);
		if (is_write) memaccs.back().set_write_through(is_writethrough);
		if (conf.record_coalescing) {
			coalescing.emplace_back();
			uint32_t lanes1 = 0;
			for (int i = 0; i < nlanes; ++i) {
				coalescing.back().set_suboff(i, offs[i]);
				if (lanes[i] > 31) std::exit(1);
				lanes1 |= ((uint32_t)1 << lanes[i]);
			}
			coalescing.back().lanes = lanes1;
		}
		if (conf.record_gwarpid) {
			gwarpids.push_back(gwarpid);
		}
// 		MI_FRAMEBUF, MI_BVH, MI_AABBS, MI_VTX_ATTRIB, MI_FACES, MI_VTX_POS;
	}
	void lmem(uint64_t off, uint8_t footprint, uint8_t sm, uint32_t lane, bool is_write, bool is_writethrough)
	{
		memaccs.push_back(MA{ (uint32_t)(off / ADDRSIZE), (uint16_t)0xffff/*, footprint*//*, sm*//*, (uint8_t)0, write_type*/ });
		memaccs.back().set_sm(sm);
		memaccs.back().set_atomic(0);
// 		memaccs.back().set_write_type(write_type);
		if (is_write) memaccs.back().set_write_through(is_writethrough);
		if (conf.record_coalescing) {
			coalescing.emplace_back();
		}
		if (conf.record_gwarpid) {
			gwarpids.push_back(0);
		}
		
	}
	void process(int nsms, int warps_per_sm, int lmem_per_thread);
	std::unordered_map<uint64_t, int> known_allocs;
	std::vector<IntAlloc> allocs;
	std::vector<int> bufidmap;
	Alloc _tmp_lmem;
	void register_allocs(const std::map<uint64_t, Alloc> &_allocs, const Alloc &lmem)
	{
		_tmp_lmem = lmem;
		allocs.clear(); known_allocs.clear(); // TODO: dont reuse allocs; we need this if we want to reuse L2, but currently this is buggy for intermediate free() calls
		bufidmap.resize(_allocs.size());
		for (const auto &a : _allocs) {
			auto i = known_allocs.insert(std::make_pair(a.first, 0));
			if (i.second) {
				std::cout << a.second.order << " is new" << std::endl;
				bufidmap[a.second.order] = i.first->second = allocs.size();
				allocs.push_back(IntAlloc{ a.first, a.second.size, total_mem(), false, a.second.ld_reqs, a.second.st_reqs, a.second.ld_txn, a.second.st_txn });
			} else {
				std::cout << a.second.order << " already in" << std::endl;
				bufidmap[a.second.order] = i.first->second;
				allocs[bufidmap[a.second.order]].size = a.second.size; 
			}
		}
		std::cout << "finished alloc reg" << std::endl;
	}
	int alloc_lmem(uint64_t size)
	{
		allocs.push_back(IntAlloc{ 0, size, total_mem(), true, _tmp_lmem.ld_reqs, _tmp_lmem.st_reqs, _tmp_lmem.ld_txn, _tmp_lmem.st_txn });
		return allocs.size() - 1;
	}
	void free_lmem()
	{
		allocs.resize(allocs.size() - 1);
	}
	uint64_t total_mem()
	{
		uint64_t ps = allocs.empty() ? 0 : allocs.back().pseudoseg + allocs.back().size;
		ps = (ps + 127) / 128 * 128;
		return ps;
	}
};


