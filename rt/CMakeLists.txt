cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(CcTrace CUDA CXX)
add_subdirectory(image)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

set(PROFILER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../profiler)

add_library(meminf SHARED meminf_default.cc)
target_include_directories(meminf PRIVATE ${PROFILER_DIR})
add_executable(rt trace.cu bvh.cc reader/ply.cc)
target_include_directories(rt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PROFILER_DIR})
target_compile_features(rt PUBLIC cxx_std_11)
target_link_libraries(rt image meminf)
target_compile_options(rt PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:SHELL:-gencode arch=compute_75,code=sm_75>")
target_link_libraries(rt image meminf)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(rt OpenMP::OpenMP_CXX)
endif()
