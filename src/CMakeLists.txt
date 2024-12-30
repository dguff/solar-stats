######################################################################
# @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
# @file        : CMakeLists
# @created     : Saturday Dec 14, 2024 21:08:05 CET
######################################################################

cmake_minimum_required(VERSION 3.12)
project(solar-stats LANGUAGES CXX VERSION 0.1)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}
  CMAKE_LIBRARY_OUTPUT_DIRECTORY)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_CONFIGURATION_TYPES Debug Release)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
if (CMAKE_BUILD_TYPE STREQUAL "Debug") 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
  message(STATUS "DEBUG: Adding '-g' option for gdb debugging")
endif()

#----------------------------------------------------------------------------
# Define few SoLAr-stats options
option(SLAR_PROFILE "Enbale profiling with gprof" OFF)

if (SLAR_PROFILE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  message(STATUS "Building solar_stats with gprof profiling feature")
endif()

#----------------------------------------------------------------------------
# Find dependencies (ROOT, RapidJSON)
find_package(ROOT REQUIRED COMPONENTS RIO Core)
find_package(RapidJSON REQUIRED)
find_package(Prob3plusplus REQUIRED)
find_library(MARLEY
  NAMES MARLEY
  PATHS "${MARLEY_DIR}/build"
  REQUIRED
  NO_DEFAULT_PATH
  )
find_library(MARLEY_ROOT
  NAMES MARLEY_ROOT
  PATHS "${MARLEY_DIR}/build"
  REQUIRED
  NO_DEFAULT_PATH
  )


include(${ROOT_USE_FILE})
include_directories(
  ${RapidJSON_INCLUDE_DIR}/..
  ${ROOT_INCLUDE_DIRS}
  ${Prob3plusplus_INCLUDE_DIR}
  ${MARLEY_DIR}/include
)

#----------------------------------------------------------------------------
# Set the lists for headers and sources
set(solar-stats-sources "")
set(solar-stats-headers "")
# Set the installation directories and paths 
set(SSTATS_CORE_DIR    "${PROJECT_SOURCE_DIR}/core")
set(SSTATS_HISTFIT_DIR "${PROJECT_SOURCE_DIR}/hist-fit")
set(SSTATS_LIB_DIR     "${CMAKE_INSTALL_PREFIX}/lib")
set(SSTATS_BIN_DIR     "${CMAKE_INSTALL_PREFIX}/bin")
set(SSTATS_INC_DEST    "${CMAKE_INSTALL_PREFIX}/include")

#----------------------------------------------------------------------------
# Compile core and hist-fit libraries
add_subdirectory( core )
add_subdirectory( hist-fit )

#----------------------------------------------------------------------------
# Installation instructions



