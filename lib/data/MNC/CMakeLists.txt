cmake_minimum_required(VERSION 3.10)
project(MNC)

find_package(PythonInterp 3.5 REQUIRED)


find_package(PythonInterp 3.5 REQUIRED)
# This goes after, since it uses PythonInterp as hint
find_package(PythonLibs 3.5 REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# This comes to hand if we also need to use the NumPy C API
exec_program(${PYTHON_EXECUTABLE}
             ARGS "-c \"import numpy; print(numpy.get_include())\""
             OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
             RETURN_VALUE NUMPY_NOT_FOUND
            )
if(NUMPY_NOT_FOUND)
    message(FATAL_ERROR "NumPy headers not found")
endif()

INCLUDE_DIRECTORIES(${NUMPY_INCLUDE_DIR})
add_subdirectory(${PROJECT_SOURCE_DIR}/pybind11 ${PROJECT_BINARY_DIR}/pybind11/bin)
add_subdirectory(${PROJECT_SOURCE_DIR}/xtl ${PROJECT_BINARY_DIR}/xtl/bin)
add_subdirectory(${PROJECT_SOURCE_DIR}/xsimd ${PROJECT_BINARY_DIR}/xsimd/bin)
add_subdirectory(${PROJECT_SOURCE_DIR}/xtensor ${PROJECT_BINARY_DIR}/xtensor/bin)
add_subdirectory(${PROJECT_SOURCE_DIR}/xtensor-python ${PROJECT_BINARY_DIR}/xtensor-python/bin)

include_directories(${PROJECT_SOURCE_DIR}/xtensor/include)
include_directories(${PROJECT_SOURCE_DIR}/xtl/include)
include_directories(${PROJECT_SOURCE_DIR}/xsimd/include)
include_directories(${PROJECT_SOURCE_DIR}/pybind11/include)
include_directories(${PROJECT_SOURCE_DIR}/xtensor-python/include)

set(PYBIND11_CPP_STANDARD -std=c++14)

pybind11_add_module(MNC MODULE MNC_py.cpp radford/mod2mat.cpp radford/mod2mat.h)

