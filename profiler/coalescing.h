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

#include <map>
#include <algorithm>

#define GETPRED(m, i) ((m >> (31 - i)) & 1)
struct Map128Elm {
	uint32_t footprint;
	uint32_t norig;
	int origs[32];
};
static const int Clinesize = 32;
static const int Csubsize = 32;
int coalesce(const uint64_t *src, int bytes_src, uint32_t preds, uint64_t *dst, uint8_t *footprints, uint32_t *off_orig, uint32_t *origidx)
{
	uint32_t tag;
	int j = 0;
	int nsubs = Clinesize / Csubsize;
	std::map<uint64_t, Map128Elm> addrmap;
	for (int i = 0; i < 32; ++i) {
		if (!GETPRED(preds, i)) continue;
		if (bytes_src > Csubsize - src[i] % Csubsize) { // unaligned
			if (bytes_src > Csubsize - src[i] % Csubsize + Csubsize)
			tag = 1 << ((src[i] / Csubsize + 1) % nsubs);
			auto ins = addrmap.insert(std::make_pair((src[i] / Csubsize + 1) / nsubs, Map128Elm{ tag, 1, { i } }));
			if (!ins.second) {
				ins.first->second.footprint |= tag;
				ins.first->second.origs[ins.first->second.norig++] = i;
			}
			continue;
		}
		tag = 1 << ((src[i] / Csubsize) % nsubs);
		auto ins = addrmap.insert(std::make_pair(src[i] / Clinesize, Map128Elm{ tag, 1, { i } }));
		if (!ins.second) {
			ins.first->second.footprint |= tag;
			ins.first->second.origs[ins.first->second.norig++] = i;
		}
	}
	off_orig[0] = 0;
	if (addrmap.size() > 32) {
		std::exit(1);
	}
	for (const auto &a : addrmap) {
		uint64_t base_addr = a.first;
		uint32_t footprint = a.second.footprint;
		dst[j] = base_addr * Clinesize;
		footprints[j] = footprint;
		off_orig[j + 1] = off_orig[j] + a.second.norig;
		std::copy(a.second.origs, a.second.origs + a.second.norig, origidx + off_orig[j]);
		++j;
	}
	return j;
} 
#undef GETPRED
