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

#include <cstdint>

extern "C" {

enum MeminfDesc {
	MI_FRAMEBUF, MI_BVH, MI_AABBS, MI_VTX_ATTRIB, MI_FACES, MI_VTX_POS, MI_SPHERES, // reserved for max' RT
	_MI_RESERVED0 = 127,
	MI_RESNET_L1_NEURONS, MI_RESNET_L1_WEIGHTS, MI_RESNET_L2_NEURONS, // some other stuff
	_MI_RESERVED1 = 4095,
	MI_VOLREN_XYZABC, // go ahead
	_MI_MAX
};
enum MeminfArg {
	MI_ARG_CAMERA, _MI_ARG_MAX
};

void meminf_name(int desc, const char *name); // this is optional, but helpful
void meminf_describe(void *ptr, int desc, int pitch = 1);
void meminf_arg(enum MeminfArg type, const char *data, int size);

}
