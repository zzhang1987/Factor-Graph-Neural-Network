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
add_subdirectory(${PROJECT_SOURCE_DIR}/3rdparty/pybind11)
add_subdirectory(${PROJECT_SOURCE_DIR}/3rdparty/xtl)
add_subdirectory(${PROJECT_SOURCE_DIR}/3rdparty/xtensor)
add_subdirectory(${PROJECT_SOURCE_DIR}/3rdparty/xtensor-python)

include_directories(${PROJECT_SOURCE_DIR}/3rdparty/xtensor/include)
include_directories(${PROJECT_SOURCE_DIR}/3rdparty/xtl/include)
include_directories(${PROJECT_SOURCE_DIR}/3rdparty/pybind11/include)
include_directories(${PROJECT_SOURCE_DIR}/3rdparty/xtensor-python/include)

set(PYBIND11_CPP_STANDARD -std=c++14)
pybind11_add_module(MNC MODULE MNC_py.cpp zb2x.h zb2x.cpp radford/mod2mat.cpp radford/mod2mat.h bnd/bnd.h bnd/bnd.cpp  ansi/cmatrix.h ansi/cmatrix.cpp ansi/nrutil.h ansi/nrutil.cpp ansi/r.cpp ansi/r.h ansi/rand2.h ansi/rand2.cpp)
