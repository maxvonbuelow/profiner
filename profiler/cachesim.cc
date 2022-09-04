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

#include "cachesim.h"

void process_parallel(MA *memaccs, uint32_t n, const IntAlloc *allocs, int lmembuf, HitCounter **rd1, HitCounter &rd2, int linesize1, int nsms);

void process(MA *memaccs, uint32_t n, const IntAlloc *allocs, int lmembuf, HitCounter &cs, int linesize, int selsm, int cachelevel)
{
	for (uint32_t i = 0; i < n; ++i) {
		MA &ma = memaccs[i];
		if (selsm != -1 && ma.get_sm() != selsm) continue;
		if (cachelevel == 1 && !ma.goes_to_l2()) continue;
		int bufid = ma.bufid;
		uint64_t addr = ma.off + allocs[bufid].pseudoseg;

		bool hit;
		if (ma.get_atomic()) {
			hit = false;
		} else {
			uint64_t line = addr / linesize;
			hit = cs(line);
		}
		if (cachelevel == 0) ma.set_l1hit(hit);
		else ma.set_l2hit(hit);
	}
}

void CacheSim::process(int nsms, int nwarps, int lmem_per_thread)
{
	uint32_t blocks_l1 = cacheconf.l1_capacity / 128;
	uint32_t assoc_l1 = blocks_l1; // fully
	uint32_t blocks_l2 = cacheconf.l2_capacity / L2LINE;
	uint32_t assoc_l2 = cacheconf.l2_assoc;
	bool clearcache = true;
	if (reuse_l1_storage.size() != nsms || clearcache) {
		reuse_l2 = CacheL2Type();
		reuse_l1_storage.resize(nsms);
		for (int i = 0; i < nsms; ++i) {
			reuse_l1_storage[i] = CacheType();
			reuse_l1_storage[i].init(blocks_l1, assoc_l1);
		}
		reuse_l2.init(blocks_l2, assoc_l2);
	}
	std::vector<HitCounter*> reuse_l1(nsms);
	uint32_t n = memaccs.size();
	uint64_t lmem = lmem_per_thread * 32 * nwarps;
	int lmembuf = -1;
	if (lmem) {
		lmembuf = alloc_lmem(lmem);
	}

	std::cout << "Local memory: " << lmem << std::endl;

	for (int i = 0; i < allocs.size(); ++i) {
		std::cout << i << " " << allocs[i].pseudoseg << "+" << allocs[i].size << std::endl;
	}
	for (uint32_t i = 0; i < n; ++i) {
		MA &ma = memaccs[i];
		int bufid;
		if (ma.bufid == 0xffff) {
			bufid = lmembuf;
		} else {
			bufid = bufidmap[ma.bufid];
		}
		ma.bufid = bufid;
		uint64_t addr = ma.off + allocs[bufid].pseudoseg; 
		if (addr >= total_mem()) {std::cerr << "sf" << addr << " " << total_mem() <<" " << bufid << "(" << ma.off << ")" << std::endl; std::exit(1);}
	}
	reuse_l2.resize(total_mem() / L2LINE);
	for (int i = 0; i < nsms; ++i) {
		reuse_l1_storage[i].resize(total_mem() / L1LINE);
		reuse_l1[i] = &reuse_l1_storage[i];
	}
	// the following is for l1
	std::cout << "Do L1 in parallel" << std::endl;
	process_parallel(memaccs.data(), n, allocs.data(), lmembuf, reuse_l1.data(), reuse_l2, L1LINE, nsms);
	// now l2
	std::cout << "Do L2" << std::endl;
	::process(memaccs.data(), n, allocs.data(), lmembuf, reuse_l2, L2LINE, -1, 1);
}

void process_parallel(MA *memaccs, uint32_t n, const IntAlloc *allocs, int lmembuf, HitCounter **rd1, HitCounter &rd2, int linesize1, int nsms)
{
	std::cout << "Process in parallel " << nsms << std::endl;
#pragma omp parallel for
	for (int i = 0; i < nsms; ++i) {
		process(memaccs, n, allocs, lmembuf, *rd1[i], linesize1, i, 0);
	}
}
