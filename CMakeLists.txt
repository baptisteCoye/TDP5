cmake_minimum_required(VERSION 2.8.11)
project(TDP2_PRCD)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -lm")

option(BUILD_TESTS "Set to ON to build the functionalities tests" ON)
option(SAVE_RESULTS "Set to ON to save the results of the program at each iteration" OFF)

set(PROJET_VERBOSE 0 CACHE INT "verbose. 0 : no information; 1 : some input and output; 2 : all data")
set(VERBOSE ON)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

add_subdirectory(src)

if(BUILD_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/src/TDPConfig.h.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/src/TDPConfig.h  )