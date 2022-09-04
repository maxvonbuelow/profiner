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

typedef __attribute__((packed)) struct {
    uint8_t sm; // WARNING: future?!
    uint16_t opcode_id; // WARNING: 16 bit maybe could cause trouble if you have extremly much memory references
    uint32_t addrs[32];
	uint64_t base;
	uint32_t preds;
	uint64_t gwarpid;
	uint64_t a(int i)
	{
		return base + addrs[i];
	}
} mem_access_t;

struct __attribute__((packed)) memop {
	enum MemorySpace { LOCAL, GENERIC, GLOBAL, SHARED, GLOBAL_TO_SHARED, SURFACE, TEXTURE, };
	char ld, st;
	MemorySpace ms;
	int s; // 1 2 4 8 16
	bool atomic;
};
