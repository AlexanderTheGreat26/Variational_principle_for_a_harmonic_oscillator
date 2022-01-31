cmake_minimum_required(VERSION 3.21)
project(Variational_principle_debug)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -ffast-math")

add_executable(Variational_principle_debug main.cpp)
