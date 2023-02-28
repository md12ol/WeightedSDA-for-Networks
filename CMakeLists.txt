cmake_minimum_required(VERSION 3.13.3)
project(SDWS)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-g")

include_directories(BitSprayer)
include_directories(Epi)

add_executable(Testy
        Epi/main.cpp
        Epi/setu.cpp
        Epi/setu.h
        Epi/stat.cpp
        Epi/filesystem.hpp
        Epi/Graph.cpp
        Epi/stat.h
        Epi/Graph.h
        BitSprayer/Bitsprayer.cpp
        BitSprayer/Bitsprayer.h Epi/WeightedGraph.h)

add_executable(SideTesty
        Epi/Graph.cpp
        Epi/Graph.h
        BitSprayer/Bitsprayer.cpp
        BitSprayer/Bitsprayer.h)