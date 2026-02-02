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

add_executable(wdeck_noise example_wdeck_noise_analysis.cc)
target_link_libraries(wdeck_noise ${ROOT_LIBRARIES} TWDeck )

install(TARGETS wdeck_smoothing  wdeck_wiener wdeck_producer wdeck_noise
  LIBRARY DESTINATION "${WDECK_BIN_DIR}"
  ARCHIVE DESTINATION "${WDECK_BIN_DIR}"
  RUNTIME DESTINATION "${WDECK_BIN_DIR}"
)

if (BUILD_TESTING)
  add_test(NAME ExampleTestSmoothing COMMAND wdeck_smoothing)
  set_tests_properties(ExampleTestSmoothing PROPERTIES 
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/examples_output_plots")

  add_test(NAME ExampleTestProducer COMMAND wdeck_producer)
  set_tests_properties(ExampleTestProducer PROPERTIES 
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/examples_output_plots")
  
  add_test(NAME ExampleTestWiener COMMAND wdeck_wiener)
  set_tests_properties(ExampleTestWiener PROPERTIES 
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/examples_output_plots")
  
  add_test(NAME ExampleTestNoise COMMAND wdeck_noise)
  set_tests_properties(ExampleTestNoise PROPERTIES 
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/examples_output_plots")
endif()
