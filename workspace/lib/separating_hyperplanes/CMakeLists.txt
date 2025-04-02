cmake_minimum_required(VERSION 3.16)
project(separating_hyperplanes)
set_property(GLOBAL PROPERTY separating_hyperplanes_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../cmake ${CMAKE_MODULE_PATH})

# add local library
get_property(is_math_ADDED GLOBAL PROPERTY math_ADDED)
if(NOT is_math_ADDED)
    add_subdirectory(../math ${CMAKE_CURRENT_BINARY_DIR}/math)
endif()

get_property(is_qpcpp_ADDED GLOBAL PROPERTY qpcpp_ADDED)
if(NOT is_qpcpp_ADDED)
    add_subdirectory(../qpcpp ${CMAKE_CURRENT_BINARY_DIR}/qpcpp)
endif()

# add libraries
add_library(separating_hyperplanes_SVM SHARED src/SVM.cpp)
target_include_directories(separating_hyperplanes_SVM PUBLIC include)
target_link_libraries(separating_hyperplanes_SVM PUBLIC
        math_Types
        qpcpp_Problem
        qpcpp_solvers_CPLEX)

add_library(separating_hyperplanes_Voronoi SHARED src/Voronoi.cpp)
target_include_directories(separating_hyperplanes_Voronoi PUBLIC include)
target_link_libraries(separating_hyperplanes_Voronoi PUBLIC
        math_Types)

# Combined library that includes both SVM and Voronoi
add_library(separating_hyperplanes INTERFACE)
target_link_libraries(separating_hyperplanes INTERFACE
        separating_hyperplanes_SVM
        separating_hyperplanes_Voronoi)

# Testing configuration
option(SEPARATING_HYPERPLANES_BUILD_TESTS "Build tests for separating_hyperplanes library" ON)

if(SEPARATING_HYPERPLANES_BUILD_TESTS)
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
        message(WARNING "Skipping tests for separating_hyperplanes library as GoogleTest is not available")
    endif()
endif()