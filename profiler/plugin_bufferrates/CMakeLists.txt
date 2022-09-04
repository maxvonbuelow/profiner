cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

project(ProfilerBufferRates)

add_library(bufferrates SHARED bufferrates.cc)
target_include_directories(bufferrates PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
