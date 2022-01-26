######################################################################
# @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
# @file        : CMakeLists
# @created     : mercoledì gen 26, 2022 11:15:32 CET
######################################################################
message(STATUS "WDECK_INCLUDE_DIR: ${WDECK_INCLUDE_DIR}")

set(HEADER_LIST 
  "${WDECK_INCLUDE_DIR}/TWDeckWfm.h"
  "${WDECK_INCLUDE_DIR}/TWDeckWfmFilter.h"
  "${WDECK_INCLUDE_DIR}/TWDeck.h"
  )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create shared libraries
ROOT_GENERATE_DICTIONARY(TWDeckDict
  ${WDECK_INCLUDE_DIR}/TWDeckWfm.h
  ${WDECK_INCLUDE_DIR}/TWDeckWfmFilter.h
  ${WDECK_INCLUDE_DIR}/TWDeck.h
  LINKDEF
  ${WDECK_INCLUDE_DIR}/TWDeckLinkDef.h
  )

add_library(TWDeck
  SHARED
  ${WaveDeck_SOURCE_DIR}/src/src/TWDeckWfm.cc 
  ${WaveDeck_SOURCE_DIR}/src/src/TWDeckWfmFilter.cc 
  ${WaveDeck_SOURCE_DIR}/src/src/TWDeck.cc
  TWDeckDict.cxx
  )

target_link_libraries(TWDeck ${ROOT_LIBRARIES})

set(TWDeck_resources 
  ${CMAKE_CURRENT_BINARY_DIR}/TWDeckDict.cxx 
  ${CMAKE_CURRENT_BINARY_DIR}/libTWDeckDict.rootmap 
  ${CMAKE_CURRENT_BINARY_DIR}/libTWDeckDict_rdict.pcm)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Install instructions
install(TARGETS TWDeck
  ARCHIVE DESTINATION        "${WDECK_LIB_DIR}"
  LIBRARY DESTINATION        "${WDECK_LIB_DIR}"
  PRIVATE_HEADER DESTINATION "${WDECK_LIB_DIR}"
  RESOURCE DESTINATION       "${WDECK_LIB_DIR}"
  RUNTIME DESTINATION        "${WDECK_BIN_DIR}")

install(FILES ${TWDeck_resources} DESTINATION "${WDECK_LIB_DIR}")

