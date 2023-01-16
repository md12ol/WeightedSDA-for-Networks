cmake_minimum_required(VERSION 3.19)
project(THADSNBit)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-g")

include_directories(BitSprayer)
include_directories(Epi)


add_executable(testarooni
        Epi/main.cpp
        Epi/setu.cpp
        Epi/setu.h
        Epi/stat.cpp
        Epi/filesystem.hpp
        Epi/stat.h )

add_executable(biterooni
        Epi/main.cpp
        Epi/setu.cpp
        Epi/setu.h
        Epi/stat.cpp
        Epi/filesystem.hpp
        Epi/stat.h
        BitSprayer/Bitsprayer.cpp
        BitSprayer/Bitsprayer.h )