cmake_minimum_required(VERSION 3.16)
project(model)
set_property(GLOBAL PROPERTY model_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../cmake ${CMAKE_MODULE_PATH})

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

# Testing configuration
# option(MODEL_BUILD_TESTS "Build tests for model library" ON)

# if(MODEL_BUILD_TESTS)
#     # Enable testing
#     enable_testing()

#     # Try to find and configure GoogleTest
#     find_package(GoogleTest)

#     # Include tests directory if GoogleTest is available
#     if(GOOGLETEST_CONFIGURED)
#         # Make GMOCK_INCLUDE_DIRS available to tests
#         include_directories(${GTEST_INCLUDE_DIRS} ${GMOCK_INCLUDE_DIRS})
#         add_subdirectory(tests)
#     else()
#         message(WARNING "Skipping tests for model library as GoogleTest is not available")
#     endif()
# endif()