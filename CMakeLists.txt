######################################################################
# @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
# @file        : CMakeLists
# @created     : mercoledì gen 26, 2022 09:27:35 CET
######################################################################

cmake_minimum_required(VERSION 3.4)
project(WaveDeck)

#------------------------------------------------------------------
# Setup TODAY
string(TIMESTAMP TODAY "%Y-%m-%d T%H:%M:%S")

#------------------------------------------------------------------
# Setup ROOT path and includes
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO) 
include(${ROOT_USE_FILE})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# configure rootlogon and project path
set(WDECK_LIB_DIR     "${PROJECT_SOURCE_DIR}/lib")
set(WDECK_BIN_DIR     "${PROJECT_SOURCE_DIR}/bin")
set(WDECK_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src/include")

#-----------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(
  ${WDECK_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR}
  )

add_definitions(${ROOT_CXX_FLAGS})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# configure rootlogon and project path
configure_file(src/include/rootlogon.C.in 
  "${PROJECT_SOURCE_DIR}/examples/rootlogon.C" @ONLY)

# configure package config
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/WaveDeckConfig.cmake.in
  cmake/WaveDeckConfig.cmake
  INSTALL_DESTINATION "${WDEC_BIN_DIR}/WaveDeck/cmake"
  PATH_VARS WDECK_LIB_DIR WDECK_BIN_DIR WDECK_INCLUDE_DIR
  )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# add subdirectory to build WaveDeck libraries
add_subdirectory(src/src)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Add executables 
SET(CMAKE_INSTALL_RPATH "${WDECK_LIB_DIR}")
add_executable(wdeck_smoothing examples/example_wdeck_smoothing.cc)
target_link_libraries(wdeck_smoothing ${ROOT_LIBRARIES} TWDeck )

install(TARGETS wdeck_smoothing 
  LIBRARY DESTINATION "${WDECK_BIN_DIR}"
  ARCHIVE DESTINATION "${WDECK_BIN_DIR}"
  RUNTIME DESTINATION "${WDECK_BIN_DIR}"
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/WaveDeckConfig.cmake 
  DESTINATION ${WDECK_BIN_DIR}/WaveDeck/cmake)

