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

#include "plugin.h" 

#include "meminf.h"

struct BufAccu {
	double sdcm_sum;
	uint64_t n;
	BufAccu() : sdcm_sum(0), n(0)
	{}
	void account(double sum, uint64_t nn)
	{
		sdcm_sum += sum;
		n += nn;
	}
};
void accumulate_hit_rates(const MA *memaccs, uint32_t n, BufAccu *accus_l1, BufAccu *accus_l2)
{
	for (uint32_t i = 0; i < n; ++i) {
		const MA &ma = memaccs[i];
		int bufid = ma.bufid;
		accus_l1[bufid].account(ma.get_l1hit(), 1);
		if (ma.goes_to_l2()) accus_l2[bufid].account(ma.get_l2hit(), 1);
	}
}
void before_sim(CacheSim &cachesim)
{
	cachesim.conf.record_coalescing = true;
	cachesim.conf.record_gwarpid = true;
}
void after_sim(CacheSim &cachesim)
{
	std::vector<BufAccu> accus_l1(cachesim.allocs.size());
	std::vector<BufAccu> accus_l2(cachesim.allocs.size());
	accumulate_hit_rates(cachesim.memaccs.data(), cachesim.memaccs.size(), accus_l1.data(), accus_l2.data());
	uint64_t cnt_l1 = 0, cnt_l2 = 0;
	for (uint32_t i = 0; i < cachesim.allocs.size(); ++i) {
		if (cachesim.allocs[i].islmem) continue;
		cnt_l1 += accus_l1[i].n;
		cnt_l2 += accus_l2[i].n;
	}
	double overall_l1 = 0;
	double overall_l2 = 0;
	const char *names[] = { "image", "bvh", "aabb", "nrm", "mesh", "geom", "SPHERES" };
	for (uint32_t i = 0; i < cachesim.allocs.size(); ++i) {
		double phit_l1 = accus_l1[i].sdcm_sum / accus_l1[i].n;
		double phit_l2 = accus_l2[i].sdcm_sum / accus_l2[i].n;
		double pbuf_l1 = (double)accus_l1[i].n / cnt_l1;
		double pbuf_l2 = (double)accus_l2[i].n / cnt_l2;
		std::cout << "Buffer " << i;
		std::string name;
		if (!cachesim.allocs[i].islmem) {
		if (accus_l1[i].n) overall_l1 += phit_l1 * pbuf_l1;
		if (accus_l2[i].n) overall_l2 += phit_l2 * pbuf_l2;
			std::cout << " @" << std::hex << cachesim.allocs[i].off << std::dec << "+" << cachesim.allocs[i].size;
			auto itmi = meminfs.find(cachesim.allocs[i].off);
			if (itmi != meminfs.end()) {
				auto itmin = meminfnames.find(itmi->second.desc);
				if (itmin != meminfnames.end()) name = itmin->second;
				else name = names[itmi->second.desc];
				std::cout << " (" << name << ")";
			}
		} else {
			name = "LMEM";
			std::cout << " LMEM";
		}
		if (meminfs.empty()) {
			name = "buf" + std::to_string(i);
		}
		std::cout << std::endl;
		if (accus_l1[i].n) std::cout << "L1: " << PERCENTIFY(phit_l1) << " % @" << PERCENTIFY(pbuf_l1) << " %" << std::endl;
		else std::cout << "L1: [buffer unused in this kernel]" << std::endl;
		if (accus_l2[i].n) std::cout << "L2: " << PERCENTIFY(phit_l2) << " % @" << PERCENTIFY(pbuf_l2) << " %" << std::endl;
		else std::cout << "L2: [buffer unused in this kernel]" << std::endl;
		double efficiency = 1. - (double)(cachesim.allocs[i].ld_txn + cachesim.allocs[i].st_txn) / ((cachesim.allocs[i].ld_reqs + cachesim.allocs[i].st_reqs) * 32.);
		if ((accus_l1[i].n || accus_l2[i].n) && !name.empty()) std::cout << "EXPORT|" << PERCENTIFY(phit_l1) << "|" << PERCENTIFY(pbuf_l1) << "|" << PERCENTIFY(phit_l2) << "|" << PERCENTIFY(pbuf_l2) << "|" << PERCENTIFY(efficiency) << "|" << name << std::endl;
	}
	std::cout << "Overall" << std::endl;
	std::cout << "FINAL L1: " << PERCENTIFY(overall_l1) << " %" << std::endl;
	std::cout << "FINAL L2: " << PERCENTIFY(overall_l2) << " %" << std::endl;
}
