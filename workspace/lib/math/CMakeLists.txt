cmake_minimum_required(VERSION 3.16)
project(math)
set_property(GLOBAL PROPERTY math_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../cmake ${CMAKE_MODULE_PATH})

# find packages
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)

# add libraries
add_library(math_Types INTERFACE)
target_include_directories(math_Types INTERFACE include ${EIGEN3_INCLUDE_DIR})

add_library(math_Combinatorics SHARED src/Combinatorics.cpp)
target_include_directories(math_Combinatorics PUBLIC include)

add_library(math_Helpers SHARED src/Helpers.cpp)
target_include_directories(math_Helpers PUBLIC include ${EIGEN3_INCLUDE_DIR})

add_library(math_collision_shapes_CollisionShape SHARED src/collision_shapes/CollisionShape.cpp)
target_include_directories(math_collision_shapes_CollisionShape PUBLIC include)
target_link_libraries(math_collision_shapes_CollisionShape PUBLIC math_Types)

add_library(math_collision_shapes_AlignedBoxCollisionShape SHARED src/collision_shapes/AlignedBoxCollisionShape.cpp)
target_include_directories(math_collision_shapes_AlignedBoxCollisionShape PUBLIC include)
target_link_libraries(math_collision_shapes_AlignedBoxCollisionShape PUBLIC
        math_collision_shapes_CollisionShape)

# Testing configuration
option(MATH_BUILD_TESTS "Build tests for math library" ON)

if(MATH_BUILD_TESTS)
    # Enable testing
    enable_testing()

    # Try to find and configure GoogleTest
    find_package(GoogleTest)

    # Include tests directory if GoogleTest is available
    if(GOOGLETEST_CONFIGURED)
        # Make GMOCK_INCLUDE_DIRS available to tests
        include_directories(${GTEST_INCLUDE_DIRS} ${GMOCK_INCLUDE_DIRS})
        add_subdirectory(tests)
    else()
        message(WARNING "Skipping tests for math library as GoogleTest is not available")
    endif()
endif()