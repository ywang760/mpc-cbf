cmake_minimum_required(VERSION 3.16)
project(math)
set_property(GLOBAL PROPERTY math_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)

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