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
#include <thread>
#include <iostream>
#include <bitset>
#include <queue>
#include <array>

#include <cstdint>
#include <iostream>
#include <cmath>

// this is a more or less dynamic variant of std::bitset; I expect that it runs at least as fast as std::bitset
template <int NMAX = 2048>
struct VariableSizeBitset { // 2048 is probably enough for all cases: 2048 * 128 = 256 bytes L1 cache; iteration for bitwise operators ignores unused chunks
	typedef uint64_t internaltype;
	internaltype data[(NMAX + sizeof(internaltype) - 1) / sizeof(internaltype)] = { 0 };
	std::size_t n, nsub;

	inline void nmaxcheck()
	{
		if (n > NMAX) {
			std::cerr << "Welcome to the future: N is too big for PLRU implementation; we did some assumptions in the code --> increase NMAX in plru.h and recompile" << std::endl;
			std::exit(1);
		}
	}

	inline VariableSizeBitset(std::size_t n) : n(n), nsub((n + sizeof(internaltype) - 1) / sizeof(internaltype))
	{
		nmaxcheck();
	}
	inline VariableSizeBitset() : n(0)
	{}
	inline void resize(std::size_t s)
	{
		n = s;
		nsub = (n + sizeof(internaltype) - 1) / sizeof(internaltype);
		nmaxcheck();
	}
	inline std::size_t size() const
	{
		return n;
	}

	inline void set(int idx, bool v)
	{
		int primary = idx / sizeof(internaltype);
		int secondary = idx % sizeof(internaltype);
		data[primary] &= ~((internaltype)1 << secondary);
		data[primary] |= (internaltype)v << secondary;
	}
	inline bool operator[](int idx)
	{
		int primary = idx / sizeof(internaltype);
		int secondary = idx % sizeof(internaltype);
		return (data[primary] >> secondary) & 1;
	}
	inline VariableSizeBitset operator~() const
	{
		VariableSizeBitset n(size());
		for (int i = 0; i < nsub; ++i) {
			n.data[i] = ~data[i];
		}
		return n;
	}
	inline VariableSizeBitset operator&(const VariableSizeBitset &oth) const
	{
		VariableSizeBitset n(size());
		for (int i = 0; i < nsub; ++i) {
			n.data[i] = data[i] & oth.data[i];
		}
		return n;
	}
	inline VariableSizeBitset operator|(const VariableSizeBitset &oth) const
	{
		VariableSizeBitset n(size());
		for (int i = 0; i < nsub; ++i) {
			n.data[i] = data[i] | oth.data[i];
		}
		return n;
	}
	inline VariableSizeBitset operator^(const VariableSizeBitset &oth) const
	{
		VariableSizeBitset n(size());
		for (int i = 0; i < nsub; ++i) {
			n.data[i] = data[i] ^ oth.data[i];
		}
		return n;
	}
	inline VariableSizeBitset &operator&=(const VariableSizeBitset &oth)
	{
		for (int i = 0; i < nsub; ++i) {
			data[i] &= oth.data[i];
		}
		return *this;
	}
	inline VariableSizeBitset &operator|=(const VariableSizeBitset &oth)
	{
		for (int i = 0; i < nsub; ++i) {
			data[i] |= oth.data[i];
		}
		return *this;
	}
	inline VariableSizeBitset &operator^=(const VariableSizeBitset &oth)
	{
		for (int i = 0; i < nsub; ++i) {
			data[i] ^= oth.data[i];
		}
		return *this;
	}
	inline bool operator==(const VariableSizeBitset &oth)
	{
		for (int i = 0; i < nsub; ++i) {
			if (data[i] != oth.data[i])
				return false;
		}
		return true;
	}
};


#define HIT 0
#define MISS 1
// coarsely adopted from https://github.com/Tony080/TreePLRU/blob/master/TreePLRU.cpp

inline uint64_t addr_get_idx(uint64_t p, uint64_t nset)
{
	return p % nset;
}
inline uint64_t addr_get_tag(uint64_t p, uint64_t nset)
{
	return p / nset;
}

struct Cell {
    bool v;
    uint64_t tag;

    inline Cell() : v(false), tag(0) {}

    inline bool isHit(uint64_t tag) {
        return v && (tag == this->tag);
    }

    inline void fetch(uint64_t address, int nset) {
        tag = addr_get_tag(address, nset);
        v = true;
    }
};

struct PLRUPermutation {
	PLRUPermutation(){}
	std::vector<VariableSizeBitset<>> mask, value, next_value;
    inline PLRUPermutation(int nway) : mask(nway), value(nway), next_value(nway)
	{
		for (int i = 0; i < nway; ++i) {
			mask[i].resize(nway - 1);
			value[i].resize(nway - 1);
			next_value[i].resize(nway - 1);
		}
		struct GenElm {
			int lvl;
			VariableSizeBitset<> mask;
			VariableSizeBitset<> value;
		};
		std::queue<GenElm> q;
		q.emplace(GenElm{ 0, nway - 1, nway - 1 });
		int i = 0;
		do {
			auto x = q.front(); q.pop();
			x.mask.set(nway - 2 - i, 1);
			auto a = x, b = x;
			a.value.set(nway - 2 - i, 0);
			b.value.set(nway - 2 - i, 1);
			q.push(a);
			q.push(b);
			++i;
		} while (!q.empty() && i < nway - 1);
		i = 0;
		while (!q.empty()) {
			auto x = q.front(); q.pop();
			VariableSizeBitset<> nv = x.mask & ~x.value;
			mask[i] = x.mask;
			value[i] = x.value;
			next_value[i] = nv;
			++i;
		}
    }
};

struct Block {
    std::vector<Cell> cell;
	VariableSizeBitset<> state;
	const PLRUPermutation *perm;

	inline void init(const PLRUPermutation *_perm, int nway)
	{
		state.resize(nway - 1);
		cell.resize(nway);
		perm = _perm;
	}

    inline uint32_t *getByTag(uint64_t tag, uint32_t *pway, int nway) {
        for (int i = 0; i < nway; ++i) {
            if (cell[i].isHit(tag)) {
                *pway = i;
                return pway;
            }
        }
        return NULL;
    }

    inline int setLRU(uint64_t address, int nset, int nway) {
        int way = 0;
        for (int i = 0; i < nway; ++i) {
            if ((state & perm->mask[i]) == perm->value[i]) {
                state ^= perm->mask[i];
                way = i;
                break;
            }
        }
        cell[way].fetch(address, nset);
		return way;
    }

    inline uint32_t *get(uint64_t address, uint32_t *pway, int nset, int nway) {
        uint32_t *d = getByTag(addr_get_tag(address, nset), pway, nway);
        if (d != NULL) {
            return d;
        }
        return d;
    }

    inline int set(uint64_t address, uint32_t *way, int nset, int nway) {
        uint32_t *p = get(address, way, nset, nway);
        if (p != NULL) {
            state &= ~perm->mask[*way];
            state |= perm->next_value[*way];
            return HIT;
        } else {
            *way = setLRU(address, nset, nway);
            return MISS;
        }
    }
};

struct CachePLRU {
	std::unique_ptr<PLRUPermutation> perm;
    std::vector<Block> block;
	int nset, nway;
	inline CachePLRU()
	{}
    inline CachePLRU(int nset, int nway)
	{
		resize(nset, nway);
	}
	inline void resize(int _nset, int _nway)
	{
		nset = _nset;
		nway = _nway;
		perm.reset(new PLRUPermutation(nway));
		block.resize(nset);
		for (int i = 0; i < nset; ++i) {
			block[i].init(perm.get(), nway);
		}
	}

    inline bool access(uint64_t address, uint32_t *way) {
        Block &b = block[addr_get_idx(address, nset)];
		int miss = b.set(address, way, nset, nway);
		return !miss;
    }
}; 
