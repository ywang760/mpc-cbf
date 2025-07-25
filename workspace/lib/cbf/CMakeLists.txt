cmake_minimum_required(VERSION 3.16)
project(cbf)
set_property(GLOBAL PROPERTY cbf_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../cmake ${CMAKE_MODULE_PATH})

# add local library
get_property(is_common_ADDED GLOBAL PROPERTY common_ADDED)
if(NOT is_common_ADDED)
    add_subdirectory(../common ${CMAKE_CURRENT_BINARY_DIR}/common)
endif()

get_property(is_model_ADDED GLOBAL PROPERTY model_ADDED)
if(NOT is_model_ADDED)
    add_subdirectory(../model ${CMAKE_CURRENT_BINARY_DIR}/model)
endif()

get_property(is_qpcpp_ADDED GLOBAL PROPERTY qpcpp_ADDED)
if(NOT is_qpcpp_ADDED)
    add_subdirectory(../qpcpp ${CMAKE_CURRENT_BINARY_DIR}/qpcpp)
endif()

get_property(is_pf_ADDED GLOBAL PROPERTY pf_ADDED)
if(NOT is_pf_ADDED)
    add_subdirectory(../particle_filter ${CMAKE_CURRENT_BINARY_DIR}/particle_filter)
endif()

if(NOT TARGET cxxopts)
    add_subdirectory(../../../third_party/cxxopts
            ${CMAKE_CURRENT_BINARY_DIR}/cxxopts)
endif()

# add libraries
add_library(cbf_detail_FovCBF SHARED src/detail/FovCBF.cpp)
target_include_directories(cbf_detail_FovCBF PUBLIC include)
target_link_libraries(cbf_detail_FovCBF PUBLIC ginac cln
        math_Helpers)

add_library(cbf_detail_ConnectivityCBF SHARED src/detail/ConnectivityCBF.cpp)
target_include_directories(cbf_detail_ConnectivityCBF PUBLIC include)
target_link_libraries(cbf_detail_ConnectivityCBF PUBLIC ginac cln
        math_Helpers
        common_Logging)

add_library(cbf_optimization_CBFQPGeneratorBase SHARED src/optimization/CBFQPGeneratorBase.cpp)
target_include_directories(cbf_optimization_CBFQPGeneratorBase PUBLIC include)
target_link_libraries(cbf_optimization_CBFQPGeneratorBase PUBLIC
        qpcpp_QPOperations
        qpcpp_Problem
        math_Helpers)

add_library(cbf_optimization_FovQPGenerator SHARED src/optimization/FovQPGenerator.cpp)
target_include_directories(cbf_optimization_FovQPGenerator PUBLIC include)
target_link_libraries(cbf_optimization_FovQPGenerator PUBLIC
        cbf_optimization_CBFQPGeneratorBase
        cbf_detail_FovCBF)

add_library(cbf_optimization_ConnectivityQPGenerator SHARED src/optimization/ConnectivityQPGenerator.cpp)
target_include_directories(cbf_optimization_ConnectivityQPGenerator PUBLIC include)
target_link_libraries(cbf_optimization_ConnectivityQPGenerator PUBLIC
        cbf_optimization_CBFQPGeneratorBase
        cbf_detail_ConnectivityCBF)

add_library(cbf_controller_FovControl SHARED src/controller/FovControl.cpp)
target_include_directories(cbf_controller_FovControl PUBLIC include)
target_link_libraries(cbf_controller_FovControl PUBLIC
        cbf_optimization_FovQPGenerator
        qpcpp_solvers_CPLEX
        model_DoubleIntegratorXYYaw)

add_library(cbf_controller_ConnectivityControl SHARED src/controller/ConnectivityControl.cpp)
target_include_directories(cbf_controller_ConnectivityControl PUBLIC include)
target_link_libraries(cbf_controller_ConnectivityControl PUBLIC
        cbf_optimization_ConnectivityQPGenerator
        qpcpp_solvers_CPLEX
        model_DoubleIntegratorXYYaw)

# add examples

add_executable(cbf_examples_CBFFovControl_example examples/fov/CBFControl_example.cpp)
target_include_directories(cbf_examples_CBFFovControl_example PUBLIC include ../../../third_party/json/single_include)
target_link_libraries(cbf_examples_CBFFovControl_example PUBLIC
        model_DoubleIntegratorXYYaw
        cbf_controller_FovControl
        cxxopts
        ParticleFilter
        PFApplications
        math_Geometry
        math_Controls
        math_Random)

add_executable(cbf_examples_CBFFormationControl_example examples/connectivity/CBFFormationControl_example.cpp)
target_include_directories(cbf_examples_CBFFormationControl_example PUBLIC include ../../../third_party/json/single_include)
target_link_libraries(cbf_examples_CBFFormationControl_example PUBLIC
        model_DoubleIntegratorXYYaw
        cbf_controller_ConnectivityControl
        cxxopts
        math_Geometry
        math_Controls
        math_Random
        spdlog::spdlog_header_only)

# Testing configuration
option(CBF_BUILD_TESTS "Build tests for cbf library" ON)

if(CBF_BUILD_TESTS)
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
        message(WARNING "Skipping tests for cbf library as GoogleTest is not available")
    endif()
endif()