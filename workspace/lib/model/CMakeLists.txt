cmake_minimum_required(VERSION 3.16)
project(model)
set_property(GLOBAL PROPERTY model_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)

# add local dependency
get_property(is_math_ADDED GLOBAL PROPERTY math_ADDED)
if(NOT is_math_ADDED)
    add_subdirectory(../math ${CMAKE_CURRENT_BINARY_DIR}/math)
endif()

# add libraries
add_library(model_DoubleIntegrator SHARED src/DoubleIntegrator.cpp)
target_include_directories(model_DoubleIntegrator PUBLIC include)
target_link_libraries(model_DoubleIntegrator PUBLIC math_Types)

add_library(model_DoubleIntegratorXYYaw SHARED src/DoubleIntegratorXYYaw.cpp)
target_include_directories(model_DoubleIntegratorXYYaw PUBLIC include)
target_link_libraries(model_DoubleIntegratorXYYaw PUBLIC model_DoubleIntegrator)

# add executable
add_executable(model_DoubleIntegratorXYYaw_example examples/DoubleIntegratorXYYaw_example.cpp)
target_include_directories(model_DoubleIntegratorXYYaw_example PUBLIC include ../../../third_party/json/single_include)
target_link_libraries(model_DoubleIntegratorXYYaw_example PUBLIC model_DoubleIntegratorXYYaw)

add_executable(model_tests_model_test tests/model_test.cpp)
target_include_directories(model_tests_model_test PUBLIC include)
target_link_libraries(model_tests_model_test PUBLIC math_Types)