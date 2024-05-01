cmake_minimum_required(VERSION 3.27)
project(proj2DA CXX)

set(CMAKE_C_STANDARD 11)

include_directories(.)
include_directories(headerFiles)
include_directories(src)

add_executable(proj2DA
        headerFiles/Graph.h
        headerFiles/MutablePriorityQueue.h
        src/Graph.cpp
        main.cpp
        headerFiles/Data.h
        src/Data.cpp
)
