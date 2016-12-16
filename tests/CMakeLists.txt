cmake_minimum_required(VERSION 2.8.11)
project(TDP2_PRCD)

INCLUDE_DIRECTORIES(
  ${CMAKE_BINARY_DIR} 
  ${CMAKE_BINARY_DIR}/src 
  ${CMAKE_SOURCE_DIR} 
  ${CMAKE_SOURCE_DIR}/src 
  ${CMAKE_SOURCE_DIR}/tests/data
  )

add_executable(test_seq test_sequentiel.c)
target_link_libraries(test_seq util perf)
