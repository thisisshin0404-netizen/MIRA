cmake_minimum_required(VERSION 3.15)
project(MIRA LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(pybind11 REQUIRED CONFIG)

pybind11_add_module(MIRA src/module.cpp)