# Test executables for qpcpp library components

# Configure tests for qpcpp library
set(TEST_SOURCES
    ProblemTest.cpp
    QPOperationsTest.cpp
    # SolverTest.cpp
    CPLEXTest.cpp
)

add_executable(qpcpp_tests ${TEST_SOURCES})

target_link_libraries(qpcpp_tests
    PUBLIC
    qpcpp_Problem
    qpcpp_QPOperations
    qpcpp_solvers_Solver
    qpcpp_solvers_CPLEX
    gtest
    gtest_main
    gmock
    gmock_main
)

# Add include directories
target_include_directories(qpcpp_tests PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${GTEST_INCLUDE_DIRS}
    ${GMOCK_INCLUDE_DIRS}
)

gtest_discover_tests(qpcpp_tests)
