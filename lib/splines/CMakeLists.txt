cmake_minimum_required(VERSION 3.16)
project(splines)
set_property(GLOBAL PROPERTY splines_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../cmake ${CMAKE_MODULE_PATH})

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# add libs
get_property(is_math_ADDED GLOBAL PROPERTY math_ADDED)
if (NOT is_math_ADDED)
    add_subdirectory(../math ${CMAKE_CURRENT_BINARY_DIR}/math)
endif ()

get_property(is_qpcpp_ADDED GLOBAL PROPERTY qpcpp_ADDED)
if (NOT is_qpcpp_ADDED)
    add_subdirectory(../qpcpp ${CMAKE_CURRENT_BINARY_DIR}/qpcpp)
endif ()

# detail/BezierOperations
add_library(splines_detail_BezierOperations SHARED src/detail/BezierOperations.cpp)
target_include_directories(splines_detail_BezierOperations PUBLIC include ${EIGEN3_INCLUDE_DIR})
target_link_libraries(splines_detail_BezierOperations PUBLIC
        math_Types
        math_Combinatorics)

#get_property(is_optimizer_ADDED GLOBAL PROPERTY optimizer_ADDED)
#if(NOT is_optimizer_ADDED)
#    add_subdirectory(../optimizer ${CMAKE_CURRENT_BINARY_DIR}/optimizer)
#endif()


# curves/Bezier
add_library(splines_curves_Bezier SHARED src/curves/Bezier.cpp)
target_include_directories(splines_curves_Bezier PUBLIC include)
target_link_libraries(splines_curves_Bezier PUBLIC
        math_Types
        splines_detail_BezierOperations)

# curves/SingleParameterPiecewiseCurve
add_library(splines_curves_SingleParameterPiecewiseCurve SHARED src/curves/SingleParameterPiecewiseCurve.cpp)
target_include_directories(splines_curves_SingleParameterPiecewiseCurve PUBLIC include)
target_link_libraries(splines_curves_SingleParameterPiecewiseCurve PUBLIC
        math_Types)

# curves/Helpers
add_library(splines_curves_Helpers SHARED src/curves/Helpers.cpp)
target_include_directories(splines_curves_Helpers PUBLIC
        include
        ../../third_party/json/single_include)
target_link_libraries(splines_curves_Helpers PUBLIC
        math_Types)

# optimization/SingleParameterCurveQPOperations
add_library(splines_optimization_QPOperations SHARED src/optimization/QPOperations.cpp)
target_include_directories(splines_optimization_QPOperations PUBLIC include)
target_link_libraries(splines_optimization_QPOperations PUBLIC
        math_Types)

add_library(splines_optimization_SingleParameterCurveQPOperations SHARED src/optimization/SingleParameterCurveQPOperations.cpp)
target_include_directories(splines_optimization_SingleParameterCurveQPOperations PUBLIC include)
target_link_libraries(splines_optimization_SingleParameterCurveQPOperations PUBLIC
        splines_optimization_QPOperations)

# optimization/BezierQPOperations
add_library(splines_optimization_BezierQPOperations SHARED src/optimization/BezierQPOperations.cpp)
target_include_directories(splines_optimization_BezierQPOperations PUBLIC include)
target_link_libraries(splines_optimization_BezierQPOperations PUBLIC
        math_Types
        splines_curves_Bezier
        splines_optimization_SingleParameterCurveQPOperations)

# optimization/SingleParameterPiecewiseCurveQPGenerator
add_library(splines_optimization_SingleParameterPiecewiseCurveQPGenerator SHARED
        src/optimization/SingleParameterPiecewiseCurveQPGenerator.cpp)
target_include_directories(splines_optimization_SingleParameterPiecewiseCurveQPGenerator PUBLIC include)
target_link_libraries(splines_optimization_SingleParameterPiecewiseCurveQPGenerator PUBLIC
        math_Helpers
        splines_curves_SingleParameterPiecewiseCurve
        splines_optimization_SingleParameterCurveQPOperations
        qpcpp_Problem)

# -------------------------- TESTS -----------------------------------------------
# Test optimization/test_BezierOptimization
add_executable(splines_tests_BezierOptimization_test tests/optimization/BezierOptimization_test.cpp)
target_include_directories(splines_tests_BezierOptimization_test PUBLIC include)
target_link_libraries(splines_tests_BezierOptimization_test PUBLIC
        math_Types
        qpcpp_solvers_CPLEX
        splines_curves_Bezier
        splines_curves_SingleParameterPiecewiseCurve
        splines_curves_Helpers
        splines_optimization_BezierQPOperations
        splines_optimization_SingleParameterPiecewiseCurveQPGenerator
)
