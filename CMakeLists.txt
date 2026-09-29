cmake_minimum_required(VERSION 3.15...4.2)

if(DEFINED SKBUILD)
    set(MYVERSION ${SKBUILD_PROJECT_VERSION})
    set(MYPROJECTNAME ${SKBUILD_PROJECT_NAME})
else()
    set(MYVERSION "0.0.0")
    set(MYPROJECTNAME "solarprop")
endif()

project(
  ${MYPROJECTNAME}
  VERSION ${MYVERSION}
  LANGUAGES CXX)

message(STATUS "Building ${MYPROJECTNAME} ${MYVERSION}")

# Find the module development requirements (requires FindPython from 3.17 or
# scikit-build-core's built-in backport)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module Development.Embed)
find_package(pybind11 CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-Wall -Werror)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: '${CMAKE_BUILD_TYPE}'")

find_package(GSL REQUIRED)
find_package(OpenMP)

if(OpenMP_FOUND)
    add_compile_options(${OpenMP_CXX_FLAGS})
    add_link_options(${OpenMP_CXX_FLAGS})
endif()

add_subdirectory(src)

if(DEFINED SKBUILD)
  add_subdirectory(data)
  add_subdirectory(scripts)
endif()
