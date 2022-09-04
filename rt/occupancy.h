#pragma once

#include <iostream>
#include <algorithm>
#include <numeric>
#include <stdexcept>

int divup(int a, int b)
{
	return (a + b - 1) / b;
}
int ceil(int a, int b)
{
	return divup(a, b) * b;
}
int floor(int a, int b)
{
	return a / b * b;
}

enum Gran { WARP, BLOCK };
struct GPUVals {
	int SM_Version;
	int Threads_per_Warp;
	int Max_Warps_per_Multiprocessor;
	int Max_Threads_per_Multiprocessor;
	int Max_Thread_Blocks_per_Multiprocessor;
	int Registers_per_Multiprocessor;
	int Max_Registers_per_Thread_Block;
	int Shared_Memory_per_Multiprocessor;
	int Max_Shared_Memory_per_Block;
	int Register_allocation_unit_size;
	Gran Register_allocation_granularity;
	int Max_Registers_per_Thread;
	int Shared_Memory_allocation_unit_size;
	int Warp_allocation_granularity;
	int Maximum_Thread_Block_Size;
	int Runtime_Shared_Memory_Per_Block;
};

static const GPUVals gpus[] = {
	{ 20, 32, 48, 1536, 8, 49152, 49152, 32768, 32768, 64, WARP, 63, 128, 2, 1024, 0 },
	{ 21, 32, 48, 1536, 8, 49152, 49152, 32768, 32768, 64, WARP, 63, 128, 2, 1024, 0 },
	{ 30, 32, 64, 2048, 16, 49152, 49152, 65536, 65536, 256, WARP, 63, 256, 4, 1024, 0 },
	{ 32, 32, 64, 2048, 16, 49152, 49152, 65536, 65536, 256, WARP, 255, 256, 4, 1024, 0 },
	{ 35, 32, 64, 2048, 16, 49152, 49152, 65536, 65536, 256, WARP, 255, 256, 4, 1024, 0 },
	{ 37, 32, 64, 2048, 16, 114688, 49152, 131072, 65536, 256, WARP, 255, 256, 4, 1024, 0 },
	{ 50, 32, 64, 2048, 32, 65536, 49152, 65536, 65536, 256, WARP, 255, 256, 4, 1024, 0 },
	{ 52, 32, 64, 2048, 32, 98304, 49152, 65536, 32768, 256, WARP, 255, 256, 2, 1024, 0 },
	{ 53, 32, 64, 2048, 32, 65536, 49152, 65536, 32768, 256, WARP, 255, 256, 2, 1024, 0 },
	{ 60, 32, 64, 2048, 32, 65536, 49152, 65536, 65536, 256, WARP, 255, 256, 2, 1024, 0 },
	{ 61, 32, 64, 2048, 32, 98304, 49152, 65536, 65536, 256, WARP, 255, 256, 4, 1024, 0 },
	{ 62, 32, 64, 2048, 32, 65536, 49152, 65536, 65536, 256, WARP, 255, 256, 4, 1024, 0 },
	{ 70, 32, 64, 2048, 32, 98304, 98304, 65536, 65536, 256, WARP, 255, 256, 4, 1024, 0 },
	{ 75, 32, 32, 1024, 16, 65536, 65536, 65536, 65536, 256, WARP, 255, 256, 4, 1024, 0 },
	{ 80, 32, 64, 2048, 32, 167936, 167936, 65536, 65536, 256, WARP, 255, 128, 4, 1024, 1024 /* toolkit 11.1 */ },
	{ 86, 32, 48, 1536, 16, 102400, 102400, 65536, 65536, 256, WARP, 255, 128, 4, 1024, 1024 /* toolkit 11.1 */ },
};

inline GPUVals find_gpu(int sm)
{
	for (int i = 0; i < sizeof(gpus) / sizeof(GPUVals); ++i) if (gpus[i].SM_Version == sm) return gpus[i];
	throw std::runtime_error("Cannot load GPU data");
}

struct Occupancy {
	int active_blocks, active_warps;
	float occupancy;
};


inline Occupancy get_occupancy(int compute, int threads_per_block, int registers_per_thread, int smem)
{
	GPUVals gpu = find_gpu(compute);

	int warps_per_block=divup(threads_per_block, gpu.Threads_per_Warp);
	int warps_limit_persm=gpu.Max_Warps_per_Multiprocessor;
	int limit_blocks_due_to_warps=std::min(gpu.Max_Thread_Blocks_per_Multiprocessor,gpu.Max_Warps_per_Multiprocessor/warps_per_block);


	int regs_per_block=gpu.Register_allocation_granularity==BLOCK ?
		ceil(ceil(warps_per_block, gpu.Warp_allocation_granularity)*registers_per_thread*gpu.Threads_per_Warp,gpu.Register_allocation_unit_size) :
		warps_per_block;

	int regs_limit_persm=gpu.Register_allocation_granularity == BLOCK ?
		gpu.Max_Registers_per_Thread_Block :
		floor(gpu.Max_Registers_per_Thread_Block/ceil(registers_per_thread*gpu.Threads_per_Warp,gpu.Register_allocation_unit_size),gpu.Warp_allocation_granularity);

	int limit_blocks_due_to_regs=registers_per_thread>gpu.Max_Registers_per_Thread ? 0 : registers_per_thread>0 ?
		regs_limit_persm/regs_per_block*(gpu.Registers_per_Multiprocessor/gpu.Max_Registers_per_Thread_Block,1) :
		gpu.Max_Thread_Blocks_per_Multiprocessor;


	int smem_per_block=ceil(smem+gpu.Runtime_Shared_Memory_Per_Block,gpu.Shared_Memory_allocation_unit_size);

	int smem_limit_persm=gpu.Max_Shared_Memory_per_Block;

	int limit_blocks_due_to_smem=smem_per_block>gpu.Max_Shared_Memory_per_Block ? 0 : smem_per_block>0 ?
		gpu.Shared_Memory_per_Multiprocessor/smem_per_block :
		gpu.Max_Thread_Blocks_per_Multiprocessor;

	Occupancy o;
	o.active_blocks = std::min(limit_blocks_due_to_warps, std::min(limit_blocks_due_to_regs, limit_blocks_due_to_smem));
	o.active_warps = o.active_blocks*warps_per_block;
	o.occupancy = (float)o.active_warps / gpu.Max_Warps_per_Multiprocessor;
	return o;
}
