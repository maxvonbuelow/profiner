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

#include <vector>
#include <cmath>
#include <random>

#include "caches/lfu.h"
#include "caches/lru.h"
#include "caches/plru.h"
#include "caches/reuse_distance.h"
#include "caches/slru.h"

struct HitCounter {
	virtual bool operator()(uint64_t address) = 0;
};
struct PLRUCacheSim : HitCounter
{
	CachePLRU plru;
	uint32_t sets;
	uint32_t assoc;
	std::vector<uint8_t> tags;
	PLRUCacheSim(uint32_t blocks, uint32_t assoc, uint64_t totalmem)
	{
		init(blocks, assoc);
	}
	PLRUCacheSim()
	{}
	void init(uint32_t b, uint32_t a)
	{
		sets = b / a;
		assoc = a;
		plru.resize(sets, assoc);
		tags.resize(b, 0);
	}
	void resize(uint64_t totalmem)
	{}

	bool operator()(uint64_t address /* 32 byte address */)
	{
		uint64_t line = address / 4;
		int subline = address % 4;
		uint8_t footprint = 1 << subline;
		uint32_t way;
		bool hit = plru.access(line, &way);
		if (!hit) tags[way] = 0; // relevant for capacity misses; cold misses already 0
		hit = hit && !(~tags[way] & footprint);
		tags[way] |= footprint;

		return hit;
	}
};

struct NaiveRDCacheSim : HitCounter
{
	uint32_t blocks;
	uint32_t assoc;
	std::vector<ReuseDistance> rd;
	NaiveRDCacheSim(uint32_t blocks, uint32_t assoc)
	{
		init(blocks, assoc);
	}
	NaiveRDCacheSim()
	{}
	void init(uint32_t b, uint32_t a)
	{
		blocks = b;
		assoc = a;
		rd.resize(b / a);
	}
	void resize(uint32_t x)
	{}
	bool operator()(uint64_t address /* 64 byte address */)
	{
		return rd[address % rd.size()](address / rd.size()) < assoc;
	}
};



// WARNING: The following caches are currently NOT USED! Please review the code before using them.




struct LRUCacheSim : HitCounter // NOTE: LRU is pointless; use reuse distances instead; I keep it here for further experiments
{
	std::vector<CacheLRU> lru;
	std::vector<uint8_t> tags;
	uint32_t sets;
	uint32_t assoc;
	uint64_t tmem;
	LRUCacheSim(uint32_t blocks, uint32_t assoc, uint64_t totalmem)
	{
		init(blocks, assoc);
		resize(totalmem);
	}
	LRUCacheSim()
	{
	}
	void init(uint32_t b, uint32_t a)
	{
		b/=4;
		sets = b / a;
		assoc = a;
		lru.resize(sets);
		tags.resize(b, 0);
		for (uint32_t i = 0; i < sets; ++i) {
			lru[i].setSize(assoc, 0);
		}
	}
	void resize(uint64_t totalmem)
	{
		totalmem/=4;
	
		tmem = (totalmem + sets - 1) / sets;
		for (uint32_t i = 0; i < sets; ++i) {
			lru[i].setDim(tmem);
		}
	}

	bool operator()(uint64_t address /* 32 byte address */)
	{
		uint64_t line = address / 4;
		int subline = address % 4;
		uint8_t footprint = 1 << subline;

		uint64_t set = line % sets;
		uint64_t setline = line / sets;

// 		return 0;
		int cacheidx;
		if (setline >= tmem) std::cerr << "ERR " << setline << " " << tmem << std::endl;
		bool hit = lru[set].find(cacheidx, setline);
		int tagidx = cacheidx + set * assoc;
		if (!hit) tags[tagidx] = 0; // relevant for capacity misses; cold misses already 0
		hit = hit && !(~tags[tagidx] & footprint);
		tags[tagidx] |= /*1 << subline*/footprint;
		return hit;
// 		sum += hit;
// 		++n;
	}
};



struct PLRUL2CacheSim : HitCounter
{
	CachePLRU plru;
	uint32_t sets;
	uint32_t assoc;
	std::vector<uint8_t> tags;
	PLRUL2CacheSim(uint32_t blocks, uint32_t assoc, uint64_t totalmem)
	{
		assoc = 11 /* 5.5, see 2 on RHS */ * 1024 * 1024 / 128 / 32 / 2;
		uint32_t sets = 32;
		init(assoc * sets, assoc);
	}
	PLRUL2CacheSim()
	{
		assoc = 11 /* 5.5, see 2 on RHS */ * 1024 * 1024 / 128 / 32 / 2;
		uint32_t sets = 32;
		init(assoc * sets, assoc);
	}
	void init(uint32_t b, uint32_t a)
	{
		sets = b / a;
		assoc = a;
		plru.resize(sets, assoc);
		tags.resize(b, 0);
	}
	void resize(uint64_t totalmem)
	{}

	bool operator()(/*double &sum, uint64_t &n, */uint64_t address /* 64 byte address */)
	{
		uint64_t line = address / 4;
		int subline = address % 4;
		uint8_t footprint = 1 << subline;

		uint32_t way;
		bool hit = plru.access(address, &way);
		if (!hit) tags[way] = 0; // relevant for capacity misses; cold misses already 0
		hit = hit && !(~tags[way] & footprint);
		tags[way] |= /*1 << subline*/footprint;

		return hit;
// 		sum += hit;
// 		++n;
	}
};

struct LRUL2CacheSim : HitCounter
{
	std::vector<CacheLRU> lru;
	std::vector<uint8_t> tags;
	uint32_t sets;
	uint32_t assoc;
	uint64_t tmem;
	LRUL2CacheSim(uint32_t blocks, uint32_t assoc, uint64_t totalmem)
	{
		init(blocks, assoc);
		resize(totalmem);
	}
	LRUL2CacheSim()
	{
		uint32_t assoc = 11 /* 5.5, see 2 on RHS */ * 1024 * 1024 / 64 / 16 / 2;
		uint32_t sets = 16;
// 		uint32_t blocks = assoc * sets;
		init(assoc, sets);
	}
	void init(uint32_t a, uint32_t s)
	{
// 		std::cout << "INIT " << b << " " << a << std::endl;
// 		b/=4;
		sets = s;
		assoc = a;
		uint64_t b = s * a;
		lru.resize(sets);
		tags.resize(b, 0);
		for (uint32_t i = 0; i < sets; ++i) {
			lru[i].setSize(assoc, 0);
		}
	}
	void resize(uint64_t totalmem)
	{
// 		totalmem/=4;
// 		std::cout << sets << std::endl;
	
		tmem = (totalmem + sets - 1) / sets;
		for (uint32_t i = 0; i < sets; ++i) {
			lru[i].setDim(tmem);
		}
	}

	bool operator()(/*double &sum, uint64_t &n, */uint64_t address /* 32 byte address */)
	{
		uint64_t line = address/* / 4*/;
// 		int subline = address % 4;
// 		footprint = 1 << subline;

		uint64_t set = line % sets;
		uint64_t setline = line / sets;

// 		return 0;
		int cacheidx;
		if (setline >= tmem) std::cerr << "ERR " << setline << " " << tmem << std::endl;
		return lru[set].find(cacheidx, setline);
	}
};

#define RN 456
static std::random_device rd;
struct RandomEvictionCacheSim : HitCounter {
	struct Cell {
		bool v;
		uint64_t tag;
		inline Cell() : v(false), tag(0) {}
	};
	Cell cells[RN];

	uint32_t sets;
	uint32_t assoc;
	std::vector<uint8_t> tags;
    std::mt19937 g;
	std::uniform_int_distribution<std::mt19937::result_type> dist;
	RandomEvictionCacheSim(uint32_t blocks, uint32_t assoc, uint64_t totalmem) : RandomEvictionCacheSim()
	{
		init(blocks, assoc);
	}
	RandomEvictionCacheSim() : g(rd()), dist(0, RN - 1)
	{}
	void init(uint32_t b, uint32_t a)
	{
		b/=4;
		sets = b / a;
		assoc = a;
		tags.resize(b, 0);
	}
	void resize(uint64_t totalmem)
	{}
	bool access(uint64_t tag, uint32_t *way)
	{
		int evict = -1;
		for (int i = 0; i < RN; ++i) {
			if (cells[i].v && cells[i].tag == tag) {
				*way = i;
				return true;
			}
			if (evict == -1 && !cells[i].v) {
				evict = i;
			}
		}
		if (evict == -1) {
			evict = dist(g);
		}
		cells[evict].v = true;
		cells[evict].tag = tag;
		*way = evict;
		return false;
	}

	bool operator()(uint64_t address /* 32 byte address */)
	{
		uint64_t line = address / 4;
		int subline = address % 4;
		uint8_t footprint = 1 << subline;
		uint32_t way;
// 		std::cout << "sf?" << 0 << std::endl;
		bool hit = access(line, &way);
		if (!hit) tags[way] = 0; // relevant for capacity misses; cold misses already 0
		hit = hit && !(~tags[way] & footprint);
		tags[way] |= /*1 << subline*/footprint;

		return hit;
	}
};
