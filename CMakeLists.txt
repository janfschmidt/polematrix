# cmake Version
cmake_minimum_required(VERSION 2.8)

project(palattice CXX C)
set(CMAKE_BUILD_TYPE Release)

# Project Version
set(PROG_VERSION 0.99)

# Install Path
# set(CMAKE_INSTALL_PREFIX "/usr/local")


# CXX-Standard
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()
# CXX All Warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")



# Libraries
find_library(GSL_LIBRARY gsl)
find_library(GSLCBLAS_LIBRARY gslcblas)
if(NOT GSL_LIBRARY OR NOT GSLCBLAS_LIBRARY)
  message(FATAL_ERROR "Gnu Scientific Library GSL not found!")
endif()

find_library(ARMADILLO_LIBRARY armadillo)
if(NOT ARMADILLO_LIBRARY)
  message(FATAL_ERROR "Armadillo Library not found!")
endif()

find_library(PALATTICE_LIBRARY palattice)
if(NOT PALATTICE_LIBRARY)
  message(FATAL_ERROR "palattice Library not found! Get it at github.com/janfschmidt/palattice")
endif()


find_package(Boost 1.30.0
  COMPONENTS
  program_options
  filesystem
  system
  )

find_package (Threads)



# build 'polematrix'
add_executable(polematrix
  main.cpp
  debug.cpp
  Configuration.cpp
  Simulation.cpp
  Tracking.cpp
  TrackingTask.cpp
  RadiationModel.cpp
  Trajectory.cpp
  ResStrengths.cpp
  )
SET_TARGET_PROPERTIES(polematrix
  PROPERTIES
  VERSION ${PROG_VERSION}
  )

# link 'polematrix'
target_link_libraries(polematrix
  ${ARMADILLO_LIBRARY}
  ${CMAKE_THREAD_LIBS_INIT}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${PALATTICE_LIBRARY}
  ${Z_LIBRARY}
  ${GSL_LIBRARY}
  ${GSLCBLAS_LIBRARY}
  )




# custom target: write header with current version
add_custom_target(version
  ${CMAKE_COMMAND}
  -DSRC=version.hpp.in
  -DDST=version.hpp
  -DPROG_VERSION=${PROG_VERSION}
  -P gitversion.cmake
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  DEPENDS .git/HEAD .git/index
  COMMENT "Write current version to version.hpp"
  )
add_dependencies(polematrix version)





# install
install(TARGETS polematrix
  DESTINATION bin
  )

# uninstall
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)
add_custom_target(uninstall
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/cmake_uninstall.cmake)


