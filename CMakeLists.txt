cmake_minimum_required(VERSION 3.28)
project(proj2DA C)

set(CMAKE_C_STANDARD 11)

include_directories(.)
include_directories(headerFiles)

add_executable(proj2DA
        headerFiles/Graph.h)
