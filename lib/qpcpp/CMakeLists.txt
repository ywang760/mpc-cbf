cmake_minimum_required(VERSION 3.16)
project(qpcpp)
set_property(GLOBAL PROPERTY qpcpp_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../cmake ${CMAKE_MODULE_PATH})

# Problem
add_library(qpcpp_Problem SHARED src/Problem.cpp)
target_include_directories(qpcpp_Problem PUBLIC include)
target_link_libraries(qpcpp_Problem PUBLIC )

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

