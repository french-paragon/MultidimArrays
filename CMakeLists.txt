cmake_minimum_required(VERSION 3.5)

project(multidimarrays LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(BUILD_TEST "Build the test suits for the library" OFF)

set(MULTIDIMARRAYS_NAME "MultidimArrays")
set(MULTIDIMARRAYS_EXPORTNAME "MultidimArraysTargets")

set(CMAKE_INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR})
set(CMAKE_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME})
set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/include/MultidimArrays.h
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include/MultidimArrays/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/include/MultidimIndexManipulators.h
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include/MultidimArrays/)

add_library(${MULTIDIMARRAYS_NAME} INTERFACE)
add_library(${MULTIDIMARRAYS_NAME}::${MULTIDIMARRAYS_NAME} ALIAS ${MULTIDIMARRAYS_NAME})

target_include_directories(
          ${MULTIDIMARRAYS_NAME}
          INTERFACE
          $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include/>
          $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/>
)

target_sources(${MULTIDIMARRAYS_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include/MultidimArrays.h ${CMAKE_CURRENT_SOURCE_DIR}/include/MultidimIndexManipulators.h)

if (BUILD_TEST)
    add_executable(simple_test test/simple_test.cpp)
    target_link_libraries(simple_test ${MULTIDIMARRAYS_NAME})
endif(BUILD_TEST)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/MultidimArrays.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MultidimArrays)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/MultidimIndexManipulators.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MultidimArrays)



install(TARGETS ${MULTIDIMARRAYS_NAME}
          EXPORT ${MULTIDIMARRAYS_EXPORTNAME}
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      )


install(
  EXPORT ${MULTIDIMARRAYS_EXPORTNAME}
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}/cmake/${MULTIDIMARRAYS_NAME}
  NAMESPACE ${MULTIDIMARRAYS_NAME}::
  FILE ${MULTIDIMARRAYS_NAME}Targets.cmake # ensure the name is different than the AlineConfig.cmake file
  )

configure_package_config_file(
  "cmake/Config.cmake.in"
  "${MULTIDIMARRAYS_NAME}Config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_CMAKEDIR}/cmake/${MULTIDIMARRAYS_NAME}
  PATH_VARS
    CMAKE_INSTALL_LIBDIR
  )

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${MULTIDIMARRAYS_NAME}ConfigVersion.cmake"
  VERSION "0.0"
  COMPATIBILITY AnyNewerVersion
)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${MULTIDIMARRAYS_NAME}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${MULTIDIMARRAYS_NAME}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_CMAKEDIR}/cmake/${MULTIDIMARRAYS_NAME} )
