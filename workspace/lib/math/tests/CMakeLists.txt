# Test executables for math library components

# Configure tests for math library
set(TEST_SOURCES
    CollisionShapeTest.cpp
    AlignedBoxCollisionShapeTest.cpp
    CombinatoricsTest.cpp
)

add_executable(math_tests ${TEST_SOURCES})

target_link_libraries(math_tests
    PUBLIC
    math_Types
    math_Helpers
    math_Combinatorics
    math_collision_shapes_CollisionShape
    math_collision_shapes_AlignedBoxCollisionShape
    gtest
    gtest_main
    gmock
    gmock_main
)

# Add include directories
target_include_directories(math_tests PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${GTEST_INCLUDE_DIRS}
    ${GMOCK_INCLUDE_DIRS}
)

gtest_discover_tests(math_tests)
