cmake_minimum_required(VERSION 3.16)
project(mpc)
set_property(GLOBAL PROPERTY mpc_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)

# add local library
get_property(is_splines_ADDED GLOBAL PROPERTY splines_ADDED)
if (NOT is_splines_ADDED)
    add_subdirectory(../splines ${CMAKE_CURRENT_BINARY_DIR}/splines)
endif ()

get_property(is_model_ADDED GLOBAL PROPERTY model_ADDED)
if(NOT is_model_ADDED)
    add_subdirectory(../model ${CMAKE_CURRENT_BINARY_DIR}/model)
endif()

get_property(is_qpcpp_ADDED GLOBAL PROPERTY qpcpp_ADDED)
if(NOT is_qpcpp_ADDED)
    add_subdirectory(../qpcpp ${CMAKE_CURRENT_BINARY_DIR}/qpcpp)
endif()

get_property(is_separating_hyperplanes_ADDED GLOBAL PROPERTY separating_hyperplanes_ADDED)
if(NOT is_separating_hyperplanes_ADDED)
    add_subdirectory(../separating_hyperplanes ${CMAKE_CURRENT_BINARY_DIR}/separating_hyperplanes)
endif()

# add libraries
add_library(mpc_optimization_PiecewiseBeizerMPCQPOperations SHARED src/optimization/PiecewiseBezierMPCQPOperations.cpp)
target_include_directories(mpc_optimization_PiecewiseBeizerMPCQPOperations PUBLIC include)
target_link_libraries(mpc_optimization_PiecewiseBeizerMPCQPOperations PUBLIC
        model_DoubleIntegrator
        splines_optimization_SingleParameterPiecewiseCurveQPGenerator
        splines_optimization_BezierQPOperations)

add_library(mpc_optimization_PiecewiseBezierMPCQPGenerator SHARED src/optimization/PiecewiseBezierMPCQPGenerator.cpp)
target_include_directories(mpc_optimization_PiecewiseBezierMPCQPGenerator PUBLIC include)
target_link_libraries(mpc_optimization_PiecewiseBezierMPCQPGenerator PUBLIC
        mpc_optimization_PiecewiseBeizerMPCQPOperations)

add_library(mpc_controller_BezierMPC SHARED src/controller/BezierMPC.cpp)
target_include_directories(mpc_controller_BezierMPC PUBLIC include)
target_link_libraries(mpc_controller_BezierMPC PUBLIC
        mpc_optimization_PiecewiseBezierMPCQPGenerator
        separating_hyperplanes_Voronoi
        qpcpp_Problem
        qpcpp_solvers_CPLEX)


# add examples
add_executable(mpc_examples_BezierMPCXYYaw_example examples/BezierMPCXYYaw_example.cpp)
target_include_directories(mpc_examples_BezierMPCXYYaw_example PUBLIC include ../../third_party/json/single_include)
target_link_libraries(mpc_examples_BezierMPCXYYaw_example PUBLIC
        mpc_controller_BezierMPC
        model_DoubleIntegratorXYYaw
        math_collision_shapes_AlignedBoxCollisionShape)

