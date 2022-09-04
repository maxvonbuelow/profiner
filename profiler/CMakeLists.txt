cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

enable_language(CUDA)

project(NVBITMemTrace)

find_package(OpenMP)
find_package(NVBit REQUIRED)

add_library(mem_trace_inject_funcs OBJECT inject_funcs.cu)
add_library(mem_trace MODULE mem_trace.cu cachesim.cc)
add_library(meminf_impl SHARED meminf.cc)

set_property(TARGET mem_trace PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_source_files_properties(inject_funcs.cu PROPERTIES COMPILE_FLAGS " -maxrregcount=24 -Xptxas -astoolspatch --keep-device-functions -arch=sm_35")

target_link_libraries(mem_trace_inject_funcs NVBit)
target_link_libraries(mem_trace NVBit cuda mem_trace_inject_funcs $<TARGET_OBJECTS:meminf_impl>)
target_include_directories(mem_trace PRIVATE $<TARGET_PROPERTY:meminf_impl,INCLUDE_DIRECTORIES>)
if(OpenMP_CXX_FOUND)
    target_link_libraries(mem_trace OpenMP::OpenMP_CXX)
endif()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/plugin_bufferrates)

target_compile_options(mem_trace PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-w>)
target_compile_options(mem_trace PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-w>)
