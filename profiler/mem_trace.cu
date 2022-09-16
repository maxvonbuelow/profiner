/* Copyright (c) 2019, NVIDIA CORPORATION. All rights reserved.
 * Copyright (c) 2021, Max von Buelow, GRIS, Technical University of Darmstadt
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *	notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *	notice, this list of conditions and the following disclaimer in the
 *	documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *	contributors may be used to endorse or promote products derived
 *	from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <string>
#include <map>
#include <unordered_set>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <stack>
#include <sstream>
#include <stdexcept>
#include <random>
#include <algorithm>
#include <list>
#include <numeric>
#include <random>
#include "meminf_storage.h"
#include "cachesim.h"
#include "plugin.h"
#include "coalescing.h"

/* every tool needs to include this once */
#include "nvbit_tool.h"

/* nvbit interface file */
#include "nvbit.h"

/* for channel */
#include "utils/channel.hpp"

/* contains definition of the mem_access_t structure */
#include "hiacc.h"

#include <stdlib.h>
#include <dlfcn.h>


/* Channel used to communicate from GPU to CPU receiving thread */
#define CHANNEL_SIZE (1l << 20)
static __managed__ ChannelDev channel_dev;
static ChannelHost channel_host;

/* receiving thread and its control variables */
pthread_t recv_thread;
volatile bool recv_thread_started = false;
volatile bool recv_thread_receiving = false;

/* skip flag used to avoid re-entry on the nvbit_callback when issuing
 * flush_channel kernel call */
bool skip_flag = false;

/* global control variables for this tool */
uint32_t instr_begin_interval = 0;
uint32_t instr_end_interval = UINT32_MAX;
int verbose = 0;

/* opcode to id map and reverse map  */
std::map<std::string, int> opcode_to_id_map;
std::map<int, std::string> id_to_opcode_map;
std::map<int, memop> memops;


#define GETPRED(m, i) ((m >> (31 - i)) & 1)


std::list<void*> plugins;
CacheSim cachesim;

void nvbit_at_init() {
	setenv("CUDA_MANAGED_FORCE_DEVICE_ALLOC", "1", 1);
	GET_VAR_INT(
		instr_begin_interval, "INSTR_BEGIN", 0,
		"Beginning of the instruction interval where to apply instrumentation");
	GET_VAR_INT(
		instr_end_interval, "INSTR_END", UINT32_MAX,
		"End of the instruction interval where to apply instrumentation");
	GET_VAR_INT(verbose, "TOOL_VERBOSE", 0, "Enable verbosity inside the tool");
	std::string pad(100, '-');
	printf("%s\n", pad.c_str());

	std::string pluginpaths;
	GET_VAR_STR(pluginpaths, "PLUGINS", "Plugins as .so files, comma separated");
	std::stringstream sstream(pluginpaths);
	std::string plugin;
	while (std::getline(sstream, plugin, ',')){
		std::cout << "Loading plugin " << plugin << std::endl;
		void *handle = dlopen(plugin.c_str(), RTLD_LAZY);
		plugins.push_back(handle);
	}

	GET_VAR_INT(cachesim.cacheconf.l1_capacity, "L1SIZE", 57 * 1024, "Capacity of L1 cache in bytes");
	GET_VAR_INT(cachesim.cacheconf.l2_capacity, "L2SIZE", 11/*11*/ /* 5.5, see 2 on RHS */ * 1024 * 1024 / 2, "Capacity of L2 cache in bytes");
	GET_VAR_INT(cachesim.cacheconf.l2_assoc, "L2ASSOC", 16, "Associativity of L2 cache");

}
/* Set used to avoid re-instrumenting the same functions multiple times */
std::unordered_set<CUfunction> already_instrumented;


memop parse_op(Instr &instr)
{
	memop res;
	res.ld = instr.isLoad();
	res.st = instr.isStore();
	memop::MemorySpace map[128];
	map[(int)InstrType::MemorySpace::LOCAL] = memop::LOCAL;
	map[(int)InstrType::MemorySpace::GENERIC] = memop::GENERIC;
	map[(int)InstrType::MemorySpace::GLOBAL] = memop::GLOBAL;
	map[(int)InstrType::MemorySpace::SHARED] = memop::SHARED;
	map[(int)InstrType::MemorySpace::GLOBAL_TO_SHARED] = memop::GLOBAL_TO_SHARED;
	map[(int)InstrType::MemorySpace::SURFACE] = memop::SURFACE;
	map[(int)InstrType::MemorySpace::TEXTURE] = memop::TEXTURE;
	res.ms = map[(int)instr.getMemorySpace()];   
	res.s = instr.getSize();
	const char *op = instr.getOpcode();
	res.atomic = strncmp(op, "ATOMG", 5) == 0;
	return res;
}

uint32_t nctax, nctay, nctaz;
uint32_t ntx, nty, ntz;
std::size_t taboff;
uint64_t lmem_min_a = std::numeric_limits<uint64_t>::max();
uint64_t max_a = std::numeric_limits<uint64_t>::min();
void instrument_function_if_needed(CUcontext ctx, CUfunction func) {
	/* Get related functions of the kernel (device function that can be
	 * called by the kernel) */
	std::vector<CUfunction> related_functions =
		nvbit_get_related_functions(ctx, func);

	/* add kernel itself to the related function vector */
	related_functions.push_back(func);

	/* iterate on function */
	for (auto f : related_functions) {
		/* "recording" function was instrumented, if set insertion failed
		 * we have already encountered this function */
		if (!already_instrumented.insert(f).second) {
			continue;
		}
		const std::vector<Instr *> &instrs = nvbit_get_instrs(ctx, f);
		if (verbose) {
			printf("Inspecting function %s at address 0x%lx\n",
				   nvbit_get_func_name(ctx, f), nvbit_get_func_addr(f));
		}

		uint32_t cnt = 0;
		/* iterate on all the static instructions in the function */
		for (auto instr : instrs) {
			if (cnt < instr_begin_interval || cnt >= instr_end_interval ||
				instr->getMemorySpace() == InstrType::MemorySpace::NONE ||
				instr->getMemorySpace() == InstrType::MemorySpace::CONSTANT) {
				cnt++;
				continue;
			}

			if (opcode_to_id_map.find(instr->getOpcode()) ==
				opcode_to_id_map.end()) {
				uint32_t opcode_id = opcode_to_id_map.size();
				opcode_to_id_map[instr->getOpcode()] = opcode_id;
				const char *opstr = instr->getOpcode();
				uint32_t opstrlen = std::strlen(opstr);
				id_to_opcode_map[opcode_id] = std::string(opstr);
				memop mop = parse_op(*instr);
				memops[opcode_id] = mop;
			}

			int opcode_id = opcode_to_id_map[instr->getOpcode()];
			int mref_idx = 0;
			/* iterate on the operands */
			for (int i = 0; i < instr->getNumOperands(); i++) {
				/* get the operand "i" */
				const InstrType::operand_t* op = instr->getOperand(i);

				if (op->type == InstrType::OperandType::MREF) {
					/* insert call to the instrumentation function with its
					 * arguments */
					nvbit_insert_call(instr, "instrument_mem", IPOINT_BEFORE);
					/* predicate value */
					nvbit_add_call_arg_guard_pred_val(instr);
					/* opcode id */
					nvbit_add_call_arg_const_val32(instr, opcode_id);
					/* memory reference 64 bit address */
					nvbit_add_call_arg_mref_addr64(instr, mref_idx);
					/* add pointer to channel_dev*/
					nvbit_add_call_arg_const_val64(instr,
												   (uint64_t)&channel_dev);
					mref_idx++;
				}
			}
			cnt++;
		}
//		 std::exit(1);
	}
	memop mop;
	mop.ld = mop.st = 0;
// 	os.write((const char*)&mop, sizeof(memop));

// 	std::cout << nctax<<nctay<<nctaz <<std::endl;
}

__global__ void flush_channel() {
	/* push memory access with negative cta id to communicate the kernel is
	 * completed */
	mem_access_t ma;
	ma.gwarpid = -1;
	channel_dev.push(&ma, sizeof(mem_access_t));

	/* flush channel */
	channel_dev.flush();
}

uint64_t translate_lmem_addr(uint64_t a, uint32_t gwarpid, int lane, int lmem_per_thread)
{
	uint64_t elm = a / 4;
	uint64_t off = a % 4;
	uint64_t translated = gwarpid * 32 * lmem_per_thread + (elm * 32 + lane) * 4 + off;
	return translated;
}
int lmem_static_nbytes;
struct MemAccList {
	std::list<mem_access_t> data;

	void clear(uint32_t x, uint32_t y, uint32_t z)
	{
		data.clear();
	}
	void operator()(const mem_access_t &ma)
	{
		data.push_back(ma);
	}
};



std::map<uint64_t, Alloc> allocs;
int allocidx = 0;
uint64_t allocidx2ptr[4096];
int kernelidx = 0;
void proc_memaccs(MemAccList &memaccblocks, int nwarps)
{
	cachesim.clear();

	cachesim.grid(nctax, nctay, nctaz, ntx, nty, ntz);

	for (void *handle : plugins) {
		decltype(&before_sim) cb = (decltype(&before_sim))dlsym(handle, "before_sim");
		(*cb)(cachesim);
	}

	std::list<mem_access_t> &allmemaccs = memaccblocks.data;

	std::cout << "Start proc" << std::endl;
	uint64_t local_min = std::numeric_limits<uint64_t>::max(), local_max = 0;

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);

	
	uint64_t lmem = lmem_static_nbytes * 32 * nwarps/* * deviceProp.multiProcessorCount*/;

	uint64_t num_raw = 0;
	uint64_t num_coa = 0;
	uint64_t num_raw_st = 0;
	uint64_t num_coa_st = 0;
	uint64_t num_local_raw = 0;
	uint64_t num_local_coa = 0;
	uint64_t num_local_raw_st = 0;
	uint64_t num_local_coa_st = 0;
	Alloc lmemcounter;
	for (mem_access_t &ma : allmemaccs) {
// 		uint32_t ncta = ma.cta_id_x + ma.cta_id_y * nctax + ma.cta_id_z * nctax * nctay;
		uint64_t addrs[32];
		int sizes[32];
		uint8_t footprints[32];
// 					CUdeviceptr addrs[32];
// 		nwarps = std::max(nwarps, (int)ma.warp_id);
		uint32_t origs[32], origoffs[33];
// 		if (!memops[ma.opcode_id].ld) continue; // TODO TODO REQUIRED FOR L2
		bool dead = true;
		for (int i = 0; i < 32; ++i) {
			if (!GETPRED(ma.preds, i)) continue;
			dead = false;
			break;
		}
		if (dead) continue;
		uint64_t srcaddrs[32];
		if (memops[ma.opcode_id].ms == memop::LOCAL) {
			for (int j = 0; j < memops[ma.opcode_id].s; j += 4) {
				int s = std::min(4, (int)memops[ma.opcode_id].s);
				int nn = 0; 
				for (int i = 0; i < 32; ++i) {
					if (!GETPRED(ma.preds, i)) continue;
					srcaddrs[i] = ma.a(i) - lmem_min_a;
					srcaddrs[i] = translate_lmem_addr(srcaddrs[i] + j, ma.gwarpid, i, lmem_static_nbytes);
					++nn;
				}
				int na = coalesce(srcaddrs, s, ma.preds, addrs, footprints, origoffs, origs);
				if (memops[ma.opcode_id].ld) {
					num_local_raw += 32;
					num_local_coa += na;
					++lmemcounter.ld_reqs;
					lmemcounter.ld_txn += na;
				}
				if (memops[ma.opcode_id].st) {
					num_local_raw_st += 32;
					num_local_coa_st += na;
					++lmemcounter.st_reqs;
					lmemcounter.st_txn += na;
				}
				int accs = memops[ma.opcode_id].ld + memops[ma.opcode_id].st;
				for (int k = 0; k < accs; ++k) {
					for (int i = 0; i < na; ++i) {
						uint64_t off = addrs[i];
						cachesim.lmem(off, footprints[i], ma.sm, i, k == 1, false /* write back *//*   k == 1 ? 2 : 0*/);
		// 							cachesim(off, sizes[i], cache->second.order, ma.sm);
					}
				}
			}
		} else {
			for (int i = 0; i < 32; ++i) {
				srcaddrs[i] = ma.a(i);
			}

			int na = coalesce(srcaddrs, memops[ma.opcode_id].s, ma.preds, addrs, footprints, origoffs, origs);
			if (!memops[ma.opcode_id].atomic) {
				if (memops[ma.opcode_id].ld) {
					num_raw += 32;
					num_coa += na;
				}
				if (memops[ma.opcode_id].st) {
					num_raw_st += 32;
					num_coa_st += na;
				}
			}
			int firstalloc = 0;
			bool first = true;
			for (int i = 0; i < na; ++i) {
				auto cache = allocs.empty() ? allocs.end() : std::prev(allocs.upper_bound(addrs[i]));
				if (cache == allocs.end() || addrs[i] < cache->first || addrs[i] + Clinesize > cache->first + (cache->second.size + Clinesize - 1) / Clinesize * Clinesize) {
					std::cout << "Seg fault 2 " << addrs[i] << std::endl;
					for (auto a : allocs) {
						std::cout << a.first << "+" << a.second.size << " " << a.second.order << std::endl;
					}
					std::exit(1);
				}
				if (!memops[ma.opcode_id].atomic) {
					if (memops[ma.opcode_id].ld) ++cache->second.ld_txn;
					if (memops[ma.opcode_id].st) ++cache->second.st_txn;
				}
				if (first) {
					first = false;
					if (!memops[ma.opcode_id].atomic) {
						if (memops[ma.opcode_id].ld) ++cache->second.ld_reqs;
						if (memops[ma.opcode_id].st) ++cache->second.st_reqs;
					}
					firstalloc = cache->second.order;
				} else {
					if (cache->second.order != firstalloc) {
						std::cout << "Different allocs accross warp?" << std::endl;
						std::exit(1);
					}
				}


				uint64_t off = addrs[i] - cache->first;

				uint8_t offs[32];
				int noffs = 0;
				for (int j = origoffs[i]; j < origoffs[i + 1]; ++j) {
					offs[noffs++] = srcaddrs[origs[j]] - addrs[i];
				}
				if (memops[ma.opcode_id].ld)
					cachesim(off, footprints[i], cache->second.order, ma.sm, origs + origoffs[i], offs, noffs, memops[ma.opcode_id].atomic, false, false, ma.gwarpid);
				if (memops[ma.opcode_id].st)
					cachesim(off, footprints[i], cache->second.order, ma.sm, origs + origoffs[i], offs, noffs, memops[ma.opcode_id].atomic, true, true /* write through */, ma.gwarpid);
			}
		}
	}
	allmemaccs = std::list<mem_access_t>();
// 	++nwarps;
	std::cout << "NWARPS: " << nwarps << std::endl;
	std::cout << "GLOBAL Coalesce efficiency (0 is best): " <<  (double)num_coa / num_raw << " #reqs: " << num_raw / 32 << " #transacs: " << num_coa << std::endl;
	std::cout << "GLOBAL Coalesce efficiency store (0 is best): " <<  (double)num_coa_st / num_raw_st << " #reqs: " << num_raw_st / 32 << " #transacs: " << num_coa_st << std::endl;
	std::cout << "LOCAL Coalesce efficiency (0 is best): " <<  (double)num_local_coa / num_local_raw << " #reqs: " << num_local_raw / 32 << " #transacs: " << num_local_coa << std::endl;
	std::cout << "LOCAL Coalesce efficiency store (0 is best): " <<  (double)num_local_coa_st / num_local_raw_st << " #reqs: " << num_local_raw_st / 32 << " #transacs: " << num_local_coa_st << std::endl;
	std::cout << "FINAL EXPN|" << num_raw / 32 << "|" << num_raw_st / 32 << "|" << num_coa << "|" << num_coa_st << "|" << num_local_raw / 32 << "|" << num_local_raw_st / 32 << "|" << num_local_coa << "|" << num_local_coa_st << std::endl;
	std::vector<uint32_t> allocvec(allocs.size());
	std::vector<Meminf> midescs(allocs.size(), Meminf{ _MI_MAX, -1 });
	std::unordered_map<uint64_t, uint32_t> ptr2buf;
	for (const auto &a : allocs) {
		std::cout << a.first << " " << a.second.size << std::endl;
		allocvec[a.second.order] = a.second.size;
		auto miit = meminfs.find(a.first);
		if (miit != meminfs.end()) {
			midescs[a.second.order] = miit->second;
		}
		ptr2buf[a.first] = a.second.order;
	}

	cachesim.register_allocs(allocs, lmemcounter);
	cachesim.process(deviceProp.multiProcessorCount, nwarps, lmem_static_nbytes/*, allocidx2ptr, midescs.data(), allocvec.data(), allocs.size()*/);

	for (void *handle : plugins) {
		decltype(&after_sim) cb = (decltype(&after_sim))dlsym(handle, "after_sim");
// 		(*cb)(cachesim);
	}
}


// std::list<mem_access_t> allmemaccs;
MemAccList allmemaccs;
int nwarps = -1;
void nvbit_at_cuda_event(CUcontext ctx, int is_exit, nvbit_api_cuda_t cbid,
						 const char *name, void *params, CUresult *pStatus) {
	if (skip_flag) return;

	if (cbid == API_CUDA_cuLaunchKernel_ptsz ||
		cbid == API_CUDA_cuLaunchKernel) {
		cuLaunchKernel_params *p = (cuLaunchKernel_params *)params;

		if (!is_exit) {
			int nregs;
			CUDA_SAFECALL(
				cuFuncGetAttribute(&nregs, CU_FUNC_ATTRIBUTE_NUM_REGS, p->f));

			int shmem_static_nbytes;
			CUDA_SAFECALL(
				cuFuncGetAttribute(&shmem_static_nbytes,
								   CU_FUNC_ATTRIBUTE_SHARED_SIZE_BYTES, p->f));


			instrument_function_if_needed(ctx, p->f);

			nvbit_enable_instrumented(ctx, p->f, true);

			nctax = p->gridDimX;
			nctay = p->gridDimY;
			nctaz = p->gridDimZ;
			ntx = p->blockDimX;
			nty = p->blockDimY;
			ntz = p->blockDimZ;
			CUDA_SAFECALL(
				cuFuncGetAttribute(&lmem_static_nbytes,
									CU_FUNC_ATTRIBUTE_LOCAL_SIZE_BYTES, p->f));


			uint64_t lmem_base = nvbit_get_local_mem_base_addr(ctx);
			std::cout << "LMEM: " << lmem_base << " +" << lmem_static_nbytes << " per thread" << std::endl;
// 			uint32_t _vals[] = { nctax, nctay, nctaz, p->blockDimX, p->blockDimY, p->blockDimZ, nregs, shmem_static_nbytes + p->sharedMemBytes, lmem_static_nbytes };


			printf(
				"Kernel %s - grid size %d,%d,%d - block size %d,%d,%d - nregs "
				"%d - shmem %d - cuda stream id %ld\n",
				nvbit_get_func_name(ctx, p->f), p->gridDimX, p->gridDimY,
				p->gridDimZ, p->blockDimX, p->blockDimY, p->blockDimZ, nregs,
				shmem_static_nbytes + p->sharedMemBytes, (uint64_t)p->hStream);
			allmemaccs.clear(nctax, nctay, nctaz);
			allocs.clear(); allocidx = 0;
// 			nwarps = -1;
			uint32_t nb = nctax * nctay * nctaz;
			uint32_t nt = p->blockDimZ * p->blockDimX * p->blockDimY;
			uint32_t nw = (nt + 31) / 32;
			nwarps = nb * nw;
			recv_thread_receiving = true;
			lmem_min_a = std::numeric_limits<uint64_t>::max();

		} else {
			/* make sure current kernel is completed */
			cudaDeviceSynchronize();
			assert(cudaGetLastError() == cudaSuccess);

			/* make sure we prevent re-entry on the nvbit_callback when issuing
			 * the flush_channel kernel */
			skip_flag = true;

			/* issue flush of channel so we are sure all the memory accesses
			 * have been pushed */
			flush_channel<<<1, 1>>>();
			cudaDeviceSynchronize();
			assert(cudaGetLastError() == cudaSuccess);

			/* unset the skip flag */
			skip_flag = false;

			/* wait here until the receiving thread has not finished with the
			 * current kernel */
			while (recv_thread_receiving) {
				pthread_yield();
			}

			proc_memaccs(allmemaccs, nwarps);
		}
	}
}


void *recv_thread_fun(void *) 
{
	char *recv_buffer = (char *)malloc(CHANNEL_SIZE);
	while (recv_thread_started) {
		uint32_t num_recv_bytes = 0;
		if (recv_thread_receiving &&
			(num_recv_bytes = channel_host.recv(recv_buffer, CHANNEL_SIZE)) >
				0) {
			uint32_t num_processed_bytes = 0;
			while (num_processed_bytes < num_recv_bytes) {
				mem_access_t *ma =
					(mem_access_t *)&recv_buffer[num_processed_bytes];

				/* when we get this cta_id_x it means the kernel has completed
				 */
				if (ma->gwarpid == (uint64_t)-1) {
					recv_thread_receiving = false;
					break;
				}

				if (ma->opcode_id != 0xfffffffeu && (memops[ma->opcode_id].ms == memop::LOCAL || memops[ma->opcode_id].ms == memop::GLOBAL)) {
					if (memops[ma->opcode_id].ms == memop::GLOBAL) {
						for (int i = 0; i < 32; ++i) {
							if (!GETPRED(ma->preds, i)) continue;
							auto cache = allocs.empty() ? allocs.end() : std::prev(allocs.upper_bound(ma->a(i)));
							if (cache == allocs.end() || ma->a(i) < cache->first || ma->a(i) + memops[ma->opcode_id].s > cache->first + (cache->second.size + Clinesize - 1) / Clinesize * Clinesize) {
								CUdeviceptr p, s;
								if (cuPointerGetAttribute(&p, CU_POINTER_ATTRIBUTE_RANGE_START_ADDR, ma->a(i)) == CUDA_SUCCESS && cuPointerGetAttribute(&s, CU_POINTER_ATTRIBUTE_RANGE_SIZE, ma->a(i)) == CUDA_SUCCESS) {
									cache = allocs.find(p);
									if (cache != allocs.end()) { // this wont happen due to clear()
										std::cout << "AGAIN " << p << "+" << s << " @" << cache->second.order << std::endl;
										cache->second.size = s;
									} else {
										if (cache != allocs.end()) std::cout << cache->first << std::endl;
										std::cout << "ADDING " << p << "+" << s << " @" << allocidx << std::endl;
										cache = allocs.emplace(p, Alloc{ s, allocidx++, 0/*shm*/, 0, 0 }).first;
									}
								} else {
									std::cerr << "SEG FAULT " << ma->a(i) << " " << memops[ma->opcode_id].ms << std::endl;
									std::exit(1);
								}
							} else if (cache != allocs.end()) {
							}
						}
					} else if (memops[ma->opcode_id].ms == memop::LOCAL) {
						for (int i = 0; i < 32; ++i) {
							if (!GETPRED(ma->preds, i)) continue;
							lmem_min_a = std::min(lmem_min_a, ma->a(i));
						}
					}
					allmemaccs(*ma);
				}

				num_processed_bytes += sizeof(mem_access_t);
			}
		}
	}

	free(recv_buffer);


	return NULL;
}

void nvbit_at_ctx_init(CUcontext ctx)
{
	recv_thread_started = true;
	channel_host.init(0, CHANNEL_SIZE, &channel_dev, NULL);
	pthread_create(&recv_thread, NULL, recv_thread_fun, NULL);
}

void nvbit_at_ctx_term(CUcontext ctx)
{
	if (recv_thread_started) {
		recv_thread_started = false;
		pthread_join(recv_thread, NULL);
	}
	for (void *handle : plugins) {
		dlclose(handle);
	}
}
