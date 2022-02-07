######################################################################
# @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
# @file        : CMakeLists
# @created     : lunedì gen 31, 2022 12:48:04 CET
######################################################################
message(STATUS "${CMAKE_PROJECT_NAME}: Building example executables")

SET(CMAKE_INSTALL_RPATH "${WDECK_LIB_DIR}")
add_executable(wdeck_smoothing example_wdeck_smoothing.cc)
target_link_libraries(wdeck_smoothing ${ROOT_LIBRARIES} TWDeck )

add_executable(wdeck_wiener example_wdeck_wiener.cc)
target_link_libraries(wdeck_wiener ${ROOT_LIBRARIES} TWDeck )

add_executable(wdeck_producer example_wdeck_waveproducer.cc)
target_link_libraries(wdeck_producer ${ROOT_LIBRARIES} TWDeck )

install(TARGETS wdeck_smoothing  wdeck_wiener wdeck_producer
  LIBRARY DESTINATION "${WDECK_BIN_DIR}"
  ARCHIVE DESTINATION "${WDECK_BIN_DIR}"
  RUNTIME DESTINATION "${WDECK_BIN_DIR}"
)


