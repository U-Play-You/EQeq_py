cmake_minimum_required(VERSION 3.14)
project(eqeq_pybind LANGUAGES CXX)

# 使用 C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(pybind11 REQUIRED)

set(SOURCES
        eqeq_core.cpp
        eqeq_pybind.cpp
)

pybind11_add_module(eqeq ${SOURCES})

if (MSVC)
    target_compile_options(eqeq PRIVATE /W3 /permissive-)
else()
    target_compile_options(eqeq PRIVATE -Wall -Wextra -Wno-unused-parameter -O3)
    target_link_libraries(eqeq PRIVATE m)
endif()

