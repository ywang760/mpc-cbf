cmake_minimum_required(VERSION 3.16)
project(qpcpp)
set_property(GLOBAL PROPERTY qpcpp_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../cmake ${CMAKE_MODULE_PATH})

# add local libs
get_property(is_math_ADDED GLOBAL PROPERTY math_ADDED)
if (NOT is_math_ADDED)
    add_subdirectory(../math ${CMAKE_CURRENT_BINARY_DIR}/math)
endif ()

# Problem
add_library(qpcpp_Problem SHARED src/Problem.cpp)
target_include_directories(qpcpp_Problem PUBLIC include)
target_link_libraries(qpcpp_Problem PUBLIC )

# QPOperations
add_library(qpcpp_QPOperations SHARED src/QPOperations.cpp)
target_include_directories(qpcpp_QPOperations PUBLIC include)
target_link_libraries(qpcpp_QPOperations PUBLIC math_Types)

# solvers/Solver
add_library(qpcpp_solvers_Solver SHARED src/solvers/Solver.cpp)
target_include_directories(qpcpp_solvers_Solver PUBLIC include)
target_link_libraries(qpcpp_solvers_Solver PUBLIC qpcpp_Problem)

# solvers/CPLEX
find_package(CPLEX REQUIRED)
add_library(qpcpp_solvers_CPLEX SHARED src/solvers/CPLEX.cpp)
target_include_directories(qpcpp_solvers_CPLEX PUBLIC include)
target_link_libraries(qpcpp_solvers_CPLEX PUBLIC
        ilocplex
        qpcpp_Problem
        ${CMAKE_DL_LIBS}
        qpcpp_solvers_Solver)

# Testing configuration
option(QPCPP_BUILD_TESTS "Build tests for qpcpp library" ON)

if(QPCPP_BUILD_TESTS)
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
        message(WARNING "Skipping tests for qpcpp library as GoogleTest is not available")
    endif()
endif()

