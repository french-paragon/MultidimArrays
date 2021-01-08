cmake_minimum_required(VERSION 3.5)

project(multidimarrays LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(${CMAKE_PROJECT_NAME} INTERFACE)
target_sources(${CMAKE_PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/MultidimArrays.h)

add_executable(simple_test test/simple_test.cpp)
target_link_libraries(simple_test ${CMAKE_PROJECT_NAME})

target_include_directories(${CMAKE_PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/MultidimArrays.h DESTINATION include/MultidimArrays)
