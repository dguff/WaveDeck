######################################################################
# @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
# @file        : CMakeLists
# @created     : mercoledì gen 26, 2022 11:15:32 CET
######################################################################
message(STATUS "${CMAKE_PROJECT_NAME}: Building project libraries")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create shared libraries
ROOT_GENERATE_DICTIONARY(TWDeckWfmDict
  "${WDECK_INCLUDE_DIR}/TWDeckWfm.h"
  "${WDECK_INCLUDE_DIR}/TWDeckWfmAnalysis.h"
  LINKDEF
  ${WDECK_INCLUDE_DIR}/TWDeckWfmLinkDef.h
  )
add_library(TWDeckWfm
  SHARED 
  ${WaveDeck_SOURCE_DIR}/src/src/TWDeckWfm.cc 
  ${WaveDeck_SOURCE_DIR}/src/src/TWDeckWfmAnalysis.cc 
  TWDeckWfmDict.cxx)
target_link_libraries(TWDeckWfm ${ROOT_LIBRARIES})


ROOT_GENERATE_DICTIONARY(TWDeckDict
  ${WDECK_INCLUDE_DIR}/TWDeckWfmFilter.h
  ${WDECK_INCLUDE_DIR}/TWDeckWfmModel.h
  ${WDECK_INCLUDE_DIR}/TWDeck.h
  LINKDEF
  ${WDECK_INCLUDE_DIR}/TWDeckLinkDef.h
  )
add_library(TWDeck
  SHARED
  ${WaveDeck_SOURCE_DIR}/src/src/TWDeckWfmFilter.cc 
  ${WaveDeck_SOURCE_DIR}/src/src/TWDeckWfmModel.cc 
  ${WaveDeck_SOURCE_DIR}/src/src/TWDeck.cc
  TWDeckDict.cxx
  )
target_link_libraries(TWDeck TWDeckWfm ${ROOT_LIBRARIES})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Install instructions

set(TWDeck_libs
  TWDeckWfm 
  TWDeck
  )

set(TWDeck_resources 
  ${CMAKE_CURRENT_BINARY_DIR}/TWDeckWfmDict.cxx 
  ${CMAKE_CURRENT_BINARY_DIR}/libTWDeckWfmDict.rootmap 
  ${CMAKE_CURRENT_BINARY_DIR}/libTWDeckWfmDict_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/TWDeckDict.cxx 
  ${CMAKE_CURRENT_BINARY_DIR}/libTWDeckDict.rootmap 
  ${CMAKE_CURRENT_BINARY_DIR}/libTWDeckDict_rdict.pcm
  )

foreach(lib IN LISTS TWDeck_libs)
  install(TARGETS ${lib}
    EXPORT                      WaveDeckTargets
    LIBRARY DESTINATION        "${WDECK_LIB_DIR}" 
    RUNTIME DESTINATION        "${WDECK_BIN_DIR}")
endforeach()

install(FILES ${TWDeck_resources} DESTINATION "${WDECK_LIB_DIR}")

