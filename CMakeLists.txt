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
configure_file(src/include/rootlogon.C.in "${PROJECT_SOURCE_DIR}/rootlogon.C" @ONLY)

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
#add_executable(chi2_ssm chi2_ssm.C)
#add_dependencies(chi2_ssm B16dataParser)
#target_link_libraries(chi2_ssm ${ROOT_LIBRARIES} ${B16DATA_LIB} NuSol  )

#install(TARGETS chi2_ssm 
  #LIBRARY DESTINATION "${CHI2_SSM_LIB_DIR}"
  #ARCHIVE DESTINATION "${CHI2_SSM_LIB_DIR}"
  #RUNTIME DESTINATION "${CHI2_SSM_BIN_DIR}")

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/WaveDeckConfig.cmake 
    DESTINATION ${WDECK_BIN_DIR}/WaveDeck/cmake)

