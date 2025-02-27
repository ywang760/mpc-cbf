cmake_minimum_required(VERSION 3.0.2)
project(fovmpc_controller)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
        mavros
        mavros_msgs
        rospy
        nav_msgs
)

# add local libs
get_property(is_mpc_cbf_ADDED GLOBAL PROPERTY mpc_cbf_ADDED)
if(NOT is_mpc_cbf_ADDED)
    add_subdirectory(../../lib/mpc_cbf ${CMAKE_CURRENT_BINARY_DIR}/mpc_cbf)
endif()

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if your package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES test_mypkg
#  CATKIN_DEPENDS mavros rospy
#  DEPENDS system_lib
)

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
# include
  ${catkin_INCLUDE_DIRS}
)

## Declare a C++ executable
## With catkin_make all packages are built within a single CMake context
## The recommended prefix ensures that target names across packages don't collide
add_executable(control_node src/control_node.cpp)
target_include_directories(control_node PUBLIC include ../../../third_party/json/single_include)
add_dependencies(control_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(control_node PUBLIC
        ${catkin_LIBRARIES}
        mpc_cbf_controller_BezierIMPCCBF
        model_DoubleIntegratorXYYaw
        math_collision_shapes_AlignedBoxCollisionShape)

add_executable(goal_node src/goal_node.cpp)
target_include_directories(goal_node PUBLIC include ../../../third_party/json/single_include)
add_dependencies(goal_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(goal_node PUBLIC ${catkin_LIBRARIES} math_Types)

add_executable(ids_node src/ids_node.cpp)
target_include_directories(ids_node PUBLIC include ../../../third_party/json/single_include)
add_dependencies(ids_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(ids_node PUBLIC
        ${catkin_LIBRARIES}
        mpc_cbf_controller_BezierIMPCCBF
        model_DoubleIntegratorXYYaw
        math_collision_shapes_AlignedBoxCollisionShape)