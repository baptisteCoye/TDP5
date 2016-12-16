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

add_executable(test_direct test_direct.c)
target_link_libraries(test_direct util)

add_executable(test_direct_ext test_direct_ext.c)
target_link_libraries(test_direct_ext util)

add_executable(test_P2M test_P2M.c)
target_link_libraries(test_P2M util m)
