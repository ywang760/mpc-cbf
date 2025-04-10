cmake_minimum_required(VERSION 3.10)

# Test executables for separating_hyperplanes library components

# Configure test sources
set(TEST_SOURCES
    VoronoiTest.cpp
    SVMTest.cpp
)

# Create test executable
add_executable(separating_hyperplanes_tests ${TEST_SOURCES})

# Link against libraries and test frameworks
target_link_libraries(separating_hyperplanes_tests
    PUBLIC
    separating_hyperplanes_SVM
    separating_hyperplanes_Voronoi
    gtest
    gtest_main
    gmock
    gmock_main
)

# Add include directories
target_include_directories(separating_hyperplanes_tests PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${GTEST_INCLUDE_DIRS}
    ${GMOCK_INCLUDE_DIRS}
)

# Discover tests
gtest_discover_tests(separating_hyperplanes_tests)
