######################################################################
# @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
# @file        : CMakeLists
# @created     : mercoledì gen 26, 2022 09:27:35 CET
######################################################################

cmake_minimum_required(VERSION 3.4)
project(WaveDeck VERSION 0.1.0)

#------------------------------------------------------------------
# Setup TODAY
string(TIMESTAMP TODAY "%Y-%m-%d T%H:%M:%S")

#------------------------------------------------------------------
# Setup ROOT path and includes
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO) 
include(${ROOT_USE_FILE})
include(ExternalProject)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Set project paths
set(WDECK_BIN_DIR     "${CMAKE_INSTALL_PREFIX}")
set(WDECK_LIB_DIR     "${WDECK_BIN_DIR}/lib")
set(WDECK_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src/include")
set(WDECK_SOURCE_DIR  "${PROJECT_SOURCE_DIR}/src/src")
set(WDECK_DOC_DIR     "${PROJECT_SOURCE_DIR}/docs")
set(WDECK_EXAMPLE_DIR "${PROJECT_SOURCE_DIR}/examples")
set(WDECK_INSTALL_CMAKE_DIR "${WDECK_LIB_DIR}/cmake/${PROJECT_NAME}")

#-----------------------------------------------------------------
# Include headers for this project

include_directories(
  ${WDECK_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR}
  )

add_definitions(${ROOT_CXX_FLAGS})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# configure rootlogon, doxygen configuration and project path
configure_file(${WDECK_INCLUDE_DIR}/rootlogon.C.in 
  "${WDECK_EXAMPLE_DIR}/rootlogon.C" @ONLY)
configure_file(${WDECK_INCLUDE_DIR}/TWDeckPATH.h.in
  "${WDECK_INCLUDE_DIR}/TWDeckPATH.h")

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# add subdirectory to build WaveDeck libraries
add_subdirectory(${WDECK_SOURCE_DIR})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Add executables 
add_subdirectory("${WDECK_EXAMPLE_DIR}")

# Add libraries targets to the build-tree export set
export(TARGETS TWDeckWfm TWDeck
  FILE "${PROJECT_BINARY_DIR}/WaveDeckTargets.cmake")

export(PACKAGE WaveDeck)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Generate documentation using doxygen 
find_package(Doxygen
  OPTIONAL_COMPONENTS dot mscgen dia)
if (DOXYGEN_FOUND)
  ExternalProject_Add(doxygen-awesome-css
    GIT_REPOSITORY "https://github.com/jothepro/doxygen-awesome-css.git"
    GIT_TAG "v1.5.0"
    SOURCE_DIR "${WDECK_DOC_DIR}/ext"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    )
  configure_file(${PROJECT_SOURCE_DIR}/wavedeck-doc.doxy.in 
    "${WDECK_DOC_DIR}/wavedeck-doc.doxy" @ONLY)

  add_custom_target(doc_doxygen ALL
    COMMAND ${DOXYGEN_EXECUTABLE} "${WDECK_DOC_DIR}/wavedeck-doc.doxy"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating ${CMAKE_PROJECT_NAME} doxygen documentation"
    VERBATIM)
  add_dependencies(doc_doxygen doxygen-awesome-css)
endif()

# configure package config
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/WaveDeckConfig.cmake.in
  cmake/WaveDeckConfig.cmake
  INSTALL_DESTINATION "${WDECK_INSTALL_CMAKE_DIR}"
  PATH_VARS WDECK_LIB_DIR WDECK_BIN_DIR WDECK_INCLUDE_DIR WDECK_INSTALL_CMAKE_DIR
  )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/WaveDeckConfig.cmake 
  DESTINATION ${WDECK_INSTALL_CMAKE_DIR})

# configure package config
export(EXPORT WaveDeckTargets
  FILE "${WDECK_INSTALL_CMAKE_DIR}/WaveDeckTargets.cmake"
  NAMESPACE WaveDeck::
)
