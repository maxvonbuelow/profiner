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

#include "meminf.h"
#include "meminf_storage.h"
#include <string>
#include <iostream>

std::unordered_map<uint64_t, Meminf> meminfs;
std::unordered_map<int, std::string> meminfnames;
MeminfArgElm meminf_args[_MI_ARG_MAX];

void meminf_name(int desc, const char *name)
{
	meminfnames[desc] = name;
}
void meminf_describe(void *ptr, int desc, int pitch)
{
// 	std::cout << "PITCH " << pitch << std::endl;
	meminfs[(uint64_t)ptr].desc = desc;
	meminfs[(uint64_t)ptr].pitch = pitch;
}
void meminf_arg(enum MeminfArg type, const char *data, int size)
{
	meminf_args[(int)type].resize(size);
	std::copy(data, data + size, meminf_args[(int)type].begin());
// 	meminf_args[(int)type].data = data;
// 	meminf_args[(int)type].size = size;
}
