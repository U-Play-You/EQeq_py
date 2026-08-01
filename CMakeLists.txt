cmake_minimum_required(VERSION 3.18)
project(eqeq_pybind LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Python3 3.14 EXACT COMPONENTS Interpreter Development.Module REQUIRED)

# Prefer an installed pybind11 package.  Fall back to a pinned copy so a clean
# checkout can still be configured without manually installing Python packages.
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v3.0.2
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(pybind11)
endif()

# main.cpp contains PYBIND11_MODULE(eqeq, ...), so it must be built as a
# Python extension module rather than as an executable requiring main().
pybind11_add_module(eqeq MODULE main.cpp)
