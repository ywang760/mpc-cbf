# Test executables for model library components

# Configure tests for model library
set(TEST_SOURCES
    DoubleIntegratorTest.cpp
    DoubleIntegratorXYYawTest.cpp
)

add_executable(model_tests ${TEST_SOURCES})

target_link_libraries(model_tests
    PUBLIC
    model_DoubleIntegrator
    model_DoubleIntegratorXYYaw
    math_Types
    gtest
    gtest_main
    gmock
    gmock_main
)

# Add include directories
target_include_directories(model_tests PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${GTEST_INCLUDE_DIRS}
    ${GMOCK_INCLUDE_DIRS}
)

gtest_discover_tests(model_tests)
