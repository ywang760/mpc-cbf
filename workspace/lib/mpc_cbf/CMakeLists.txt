cmake_minimum_required(VERSION 3.16)
project(mpc_cbf)
set_property(GLOBAL PROPERTY mpc_cbf_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)

# add local library
get_property(is_common_ADDED GLOBAL PROPERTY common_ADDED)
if(NOT is_common_ADDED)
    add_subdirectory(../common ${CMAKE_CURRENT_BINARY_DIR}/common)
endif()

get_property(is_mpc_ADDED GLOBAL PROPERTY mpc_ADDED)
if(NOT is_mpc_ADDED)
    add_subdirectory(../mpc ${CMAKE_CURRENT_BINARY_DIR}/mpc)
endif()

get_property(is_cbf_ADDED GLOBAL PROPERTY cbf_ADDED)
if(NOT is_cbf_ADDED)
    add_subdirectory(../cbf ${CMAKE_CURRENT_BINARY_DIR}/cbf)
endif()

if(NOT TARGET cxxopts)
    add_subdirectory(../../../third_party/cxxopts
            ${CMAKE_CURRENT_BINARY_DIR}/cxxopts)
endif()

get_property(is_pf_ADDED GLOBAL PROPERTY pf_ADDED)
if(NOT is_pf_ADDED)
    add_subdirectory(../particle_filter ${CMAKE_CURRENT_BINARY_DIR}/particle_filter)
endif()

# add libraries

# optimization base classes
add_library(mpc_cbf_optimization_MPCCBFQPOperationsBase SHARED src/optimization/MPCCBFQPOperationsBase.cpp)
target_include_directories(mpc_cbf_optimization_MPCCBFQPOperationsBase PUBLIC include)
target_link_libraries(mpc_cbf_optimization_MPCCBFQPOperationsBase PUBLIC
        mpc_optimization_PiecewiseBeizerMPCQPOperations)

add_library(mpc_cbf_optimization_MPCCBFQPGeneratorBase SHARED src/optimization/MPCCBFQPGeneratorBase.cpp)
target_include_directories(mpc_cbf_optimization_MPCCBFQPGeneratorBase PUBLIC include)
target_link_libraries(mpc_cbf_optimization_MPCCBFQPGeneratorBase PUBLIC
        mpc_cbf_optimization_MPCCBFQPOperationsBase
        mpc_optimization_PiecewiseBezierMPCQPGenerator)

# FoV-specific implementation
add_library(mpc_cbf_optimization_FovMPCCBFQPOperations SHARED src/optimization/FovMPCCBFQPOperations.cpp)
target_include_directories(mpc_cbf_optimization_FovMPCCBFQPOperations PUBLIC include)
target_link_libraries(mpc_cbf_optimization_FovMPCCBFQPOperations PUBLIC
        mpc_cbf_optimization_MPCCBFQPOperationsBase
        mpc_optimization_PiecewiseBeizerMPCQPOperations
        cbf_detail_FovCBF)

add_library(mpc_cbf_optimization_FovMPCCBFQPGenerator SHARED src/optimization/FovMPCCBFQPGenerator.cpp)
target_include_directories(mpc_cbf_optimization_FovMPCCBFQPGenerator PUBLIC include)
target_link_libraries(mpc_cbf_optimization_FovMPCCBFQPGenerator PUBLIC
        mpc_cbf_optimization_MPCCBFQPGeneratorBase
        mpc_cbf_optimization_FovMPCCBFQPOperations
        mpc_optimization_PiecewiseBeizerMPCQPOperations
        mpc_optimization_PiecewiseBezierMPCQPGenerator
        cbf_detail_FovCBF)

# 本地版本的 ConnectivityCBF 实现
add_library(mpc_cbf_detail_ConnectivityCBF SHARED src/detail/ConnectivityCBF.cpp)
target_include_directories(mpc_cbf_detail_ConnectivityCBF PUBLIC include)
target_link_libraries(mpc_cbf_detail_ConnectivityCBF PUBLIC
    ginac
    cln
    math_Helpers
    common_Logging
)

# Connectivity-specific implementation
add_library(mpc_cbf_optimization_ConnectivityMPCCBFQPOperations SHARED src/optimization/ConnectivityMPCCBFQPOperations.cpp)
target_include_directories(mpc_cbf_optimization_ConnectivityMPCCBFQPOperations PUBLIC include)
target_link_libraries(mpc_cbf_optimization_ConnectivityMPCCBFQPOperations PUBLIC
        mpc_cbf_optimization_MPCCBFQPOperationsBase
        mpc_optimization_PiecewiseBeizerMPCQPOperations
        mpc_cbf_detail_ConnectivityCBF)

add_library(mpc_cbf_optimization_ConnectivityMPCCBFQPGenerator SHARED src/optimization/ConnectivityMPCCBFQPGenerator.cpp)
target_include_directories(mpc_cbf_optimization_ConnectivityMPCCBFQPGenerator PUBLIC include)
target_link_libraries(mpc_cbf_optimization_ConnectivityMPCCBFQPGenerator PUBLIC
        mpc_cbf_optimization_MPCCBFQPGeneratorBase
        mpc_cbf_optimization_ConnectivityMPCCBFQPOperations
        mpc_optimization_PiecewiseBeizerMPCQPOperations
        mpc_optimization_PiecewiseBezierMPCQPGenerator
        mpc_cbf_detail_ConnectivityCBF)

# Controllers
# MPC CBF controller (commented out - needs refactoring)
# add_library(mpc_cbf_controller_BezierMPCCBF SHARED src/controller/BezierMPCCBF.cpp)
# target_include_directories(mpc_cbf_controller_BezierMPCCBF PUBLIC include)
# target_link_libraries(mpc_cbf_controller_BezierMPCCBF PUBLIC
#         mpc_cbf_optimization_PiecewiseBezierMPCCBFQPGenerator
#         qpcpp_solvers_CPLEX)

add_library(mpc_cbf_controller_FovBezierIMPCCBF SHARED src/controller/FovBezierIMPCCBF.cpp)
target_include_directories(mpc_cbf_controller_FovBezierIMPCCBF PUBLIC include)
target_link_libraries(mpc_cbf_controller_FovBezierIMPCCBF PUBLIC
        mpc_cbf_optimization_FovMPCCBFQPGenerator
        qpcpp_solvers_CPLEX
        separating_hyperplanes_Voronoi)

add_library(mpc_cbf_controller_ConnectivityIMPCCBF SHARED src/controller/ConnectivityIMPCCBF.cpp)
target_include_directories(mpc_cbf_controller_ConnectivityIMPCCBF PUBLIC include)
target_link_libraries(mpc_cbf_controller_ConnectivityIMPCCBF PUBLIC
        mpc_cbf_optimization_ConnectivityMPCCBFQPGenerator
        qpcpp_solvers_CPLEX)

# add examples
# FoV Example using MPC CBF (commented out - needs refactoring)
# add_executable(mpc_cbf_examples_BezierMPCCBFXYYaw_example examples/fov/BezierMPCCBFXYYaw_example.cpp)
# target_include_directories(mpc_cbf_examples_BezierMPCCBFXYYaw_example PUBLIC include ../../../third_party/json/single_include)
# target_link_libraries(mpc_cbf_examples_BezierMPCCBFXYYaw_example PUBLIC
#        mpc_cbf_controller_BezierMPCCBF
#        model_DoubleIntegratorXYYaw
#        math_collision_shapes_AlignedBoxCollisionShape
#        cxxopts
#        common_Logging
#        spdlog::spdlog_header_only)

# FoV Example using Iterative MPC CBF
add_executable(mpc_cbf_examples_BezierIMPCCBFPFXYYaw_example examples/fov/BezierIMPCCBFPFXYYaw_example.cpp)
target_include_directories(mpc_cbf_examples_BezierIMPCCBFPFXYYaw_example PUBLIC include ../../../third_party/json/single_include)
target_link_libraries(mpc_cbf_examples_BezierIMPCCBFPFXYYaw_example PUBLIC
        mpc_cbf_controller_FovBezierIMPCCBF
        model_DoubleIntegratorXYYaw
        math_collision_shapes_AlignedBoxCollisionShape
        math_Geometry
        math_Random
        cxxopts
        common_Logging
        spdlog::spdlog_header_only)

# Connectivity Example using Iterative MPC CBF
add_executable(mpc_cbf_examples_MPCCBFFormationControl_example examples/connectivity/MPCCBFFormationControl_example.cpp)
target_include_directories(mpc_cbf_examples_MPCCBFFormationControl_example PUBLIC include ../../../third_party/json/single_include)
target_link_libraries(mpc_cbf_examples_MPCCBFFormationControl_example PUBLIC
        mpc_cbf_controller_ConnectivityIMPCCBF
        model_DoubleIntegratorXYYaw
        math_collision_shapes_AlignedBoxCollisionShape
        math_Geometry
        math_Random
        cxxopts
        common_Logging
        spdlog::spdlog_header_only)

target_include_directories(mpc_cbf_detail_ConnectivityCBF PUBLIC
    include
    ${CMAKE_CURRENT_SOURCE_DIR}/../math/include
)
