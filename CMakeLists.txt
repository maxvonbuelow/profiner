cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(Profiler CUDA CXX)

add_subdirectory(${CMAKE_SOURCE_DIR}/profiler)
add_subdirectory(${CMAKE_SOURCE_DIR}/rt)
