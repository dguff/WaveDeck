######################################################################
# @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
# @file        : CMakeLists.txt
# @created     : Wednesday Dec 10, 2025 12:53:40 CET
######################################################################

# Fetch GoogleTest using FetchContent (modern CMake approach)
include(FetchContent)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

# Enable GoogleTest discovery functionality
include(GoogleTest)

# Create the test executable
add_executable(TWDeckTests TWDeckUtils_test.cc TWDeckWfm_test.cc TWDeckFilter_test.cc TWDeckModel_test.cc TWDeck_test.cc)

# Link necessary libraries
target_link_libraries(TWDeckTests
  PRIVATE
  gtest_main
  TWDeckWfm 
  TWDeck
  ${ROOT_LIBRARIES}
)

# Ensure the test executable can find your headers (TWDeckWfm.h)
target_include_directories(TWDeckTests PRIVATE ${WDECK_INCLUDE_DIR})

# Register the tests with CTest so you can run them easily
gtest_discover_tests(TWDeckTests)
