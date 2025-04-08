cmake_minimum_required(VERSION 3.16)
project(cbf)
set_property(GLOBAL PROPERTY cbf_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)

# add local library
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
add_library(cbf_detail_FovCBF SHARED src/detail/cbf.cpp)
target_include_directories(cbf_detail_FovCBF PUBLIC include)
target_link_libraries(cbf_detail_FovCBF PUBLIC ginac cln
        math_Helpers)

add_library(cbf_optimization_CBFQPOperations SHARED src/optimization/CBFQPOperations.cpp)
target_include_directories(cbf_optimization_CBFQPOperations PUBLIC include)
target_link_libraries(cbf_optimization_CBFQPOperations PUBLIC
        qpcpp_QPOperations
        cbf_detail_FovCBF)

add_library(cbf_optimization_CBFQPGenerator SHARED src/optimization/CBFQPGenerator.cpp)
target_include_directories(cbf_optimization_CBFQPGenerator PUBLIC include)
target_link_libraries(cbf_optimization_CBFQPGenerator PUBLIC
        cbf_optimization_CBFQPOperations
        qpcpp_Problem
        math_Helpers)

add_library(cbf_controller_CBFControl SHARED src/controller/CBFControl.cpp)
target_include_directories(cbf_controller_CBFControl PUBLIC include)
target_link_libraries(cbf_controller_CBFControl PUBLIC
        cbf_optimization_CBFQPGenerator
        qpcpp_solvers_CPLEX
        model_DoubleIntegratorXYYaw)

# add examples
add_executable(cbf_examples_CBFControl_example examples/CBFControl_example.cpp)
target_include_directories(cbf_examples_CBFControl_example PUBLIC include ../../../third_party/json/single_include)
target_link_libraries(cbf_examples_CBFControl_example PUBLIC
        model_DoubleIntegratorXYYaw
        cbf_controller_CBFControl
        cxxopts
        ParticleFilter
        PFApplications
        math_Geometry
        math_Controls
        math_Random)

#add_executable(cbf_examples_CBFControlStaticTarget_example examples/CBFControlStaticTarget_example.cpp)
#target_include_directories(cbf_examples_CBFControlStaticTarget_example PUBLIC include ../../../third_party/json/single_include)
#target_link_libraries(cbf_examples_CBFControlStaticTarget_example PUBLIC
#        model_DoubleIntegratorXYYaw
#        cbf_controller_CBFControl)

#add_executable(cbf_examples_CBFControlSlack_example examples/CBFControlSlack_example.cpp)
#target_include_directories(cbf_examples_CBFControlSlack_example PUBLIC include ../../../third_party/json/single_include)
#target_link_libraries(cbf_examples_CBFControlSlack_example PUBLIC
#        model_DoubleIntegratorXYYaw
#        cbf_controller_CBFControl
#        ParticleFilter
#        cxxopts)
