cmake_minimum_required(VERSION 3.16)
project(particle_filter)
set_property(GLOBAL PROPERTY pf_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)

# add local dependency
get_property(is_math_ADDED GLOBAL PROPERTY math_ADDED)
if(NOT is_math_ADDED)
    add_subdirectory(../math ${CMAKE_CURRENT_BINARY_DIR}/math)
endif()

# add libraries
add_library(ParticleFilter SHARED src/detail/particle_filter.cpp)
target_include_directories(ParticleFilter PUBLIC include ${EIGEN3_INCLUDE_DIR})
# target_link_libraries(ParticleFilter PUBLIC ginac cln)

# Add the PFApplications library with implementation
add_library(PFApplications SHARED src/pf_applications.cpp)
target_include_directories(PFApplications PUBLIC include ${EIGEN3_INCLUDE_DIR})
target_link_libraries(PFApplications PUBLIC ParticleFilter math_Types math_Geometry)

# add examples
add_executable(pf_test tests/pf_test.cpp)
target_include_directories(pf_test PUBLIC include)
target_link_libraries(pf_test PUBLIC ParticleFilter)
