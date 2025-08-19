cmake_minimum_required(VERSION 3.5)

project(multidimarrays LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)

option(BUILD_TEST "Build the test suits for the library" OFF)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/include/MultidimArrays.h
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include/MultidimArrays/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/include/MultidimIndexManipulators.h
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include/MultidimArrays/)

add_library(${CMAKE_PROJECT_NAME} INTERFACE)
add_library(Multidim::Arrays ALIAS ${CMAKE_PROJECT_NAME})

target_include_directories(
          ${CMAKE_PROJECT_NAME}
          INTERFACE
          $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include/>
          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/>
)

target_sources(${CMAKE_PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include/MultidimArrays.h ${CMAKE_CURRENT_SOURCE_DIR}/include/MultidimIndexManipulators.h)

if (BUILD_TEST)
    add_executable(simple_test test/simple_test.cpp)
    target_link_libraries(simple_test ${CMAKE_PROJECT_NAME})
endif(BUILD_TEST)

target_include_directories(${CMAKE_PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/MultidimArrays.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MultidimArrays)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/MultidimIndexManipulators.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MultidimArrays)
