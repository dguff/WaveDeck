######################################################################
# @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
# @file        : CMakeLists
# @created     : mercoledì gen 26, 2022 11:15:32 CET
######################################################################
message(STATUS "${CMAKE_PROJECT_NAME}: Building project libraries")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create shared libraries
add_library(TWDeckWfm
  SHARED 
  ${WDECK_SOURCE_DIR}/TWDeckWfm.cc 
  ${WDECK_INCLUDE_DIR}/TWDeckWfm.h
  )

target_link_libraries(TWDeckWfm PUBLIC ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__TWDeckWfm
  ${WDECK_INCLUDE_DIR}/TWDeckWfm.h
  MODULE TWDeckWfm
  LINKDEF ${WDECK_INCLUDE_DIR}/TWDeckWfmLinkDef.h
  )


add_library(TWDeck
  SHARED
  ${WDECK_SOURCE_DIR}/TWDeckWfmFilter.cc 
  ${WDECK_INCLUDE_DIR}/TWDeckWfmFilter.h
  ${WDECK_SOURCE_DIR}/TWDeckWfmModel.cc 
  ${WDECK_INCLUDE_DIR}/TWDeckWfmModel.h
  ${WDECK_SOURCE_DIR}/TWDeck.cc
  ${WDECK_INCLUDE_DIR}/TWDeck.h
  )
target_link_libraries(TWDeck PUBLIC TWDeckWfm ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__TWDeck
  ${WDECK_INCLUDE_DIR}/TWDeckWfmFilter.h
  ${WDECK_INCLUDE_DIR}/TWDeckWfmModel.h
  ${WDECK_INCLUDE_DIR}/TWDeck.h
  MODULE TWDeck
  LINKDEF ${WDECK_INCLUDE_DIR}/TWDeckLinkDef.h
  )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Install instructions

set(TWDeck_libs
  TWDeckWfm 
  TWDeck
  )

set(TWDeck_resources 
  ${CMAKE_CURRENT_BINARY_DIR}/G__TWDeckWfm.cxx 
  ${CMAKE_CURRENT_BINARY_DIR}/libTWDeckWfm.rootmap 
  ${CMAKE_CURRENT_BINARY_DIR}/libTWDeckWfm_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/G__TWDeck.cxx 
  ${CMAKE_CURRENT_BINARY_DIR}/libTWDeck.rootmap 
  ${CMAKE_CURRENT_BINARY_DIR}/libTWDeck_rdict.pcm
  )

foreach(lib IN LISTS TWDeck_libs)
  install(TARGETS ${lib}
    EXPORT                      WaveDeckTargets
    LIBRARY DESTINATION        "${WDECK_LIB_DIR}" 
    RUNTIME DESTINATION        "${WDECK_BIN_DIR}")
endforeach()

install(FILES ${TWDeck_resources} DESTINATION "${WDECK_LIB_DIR}")

