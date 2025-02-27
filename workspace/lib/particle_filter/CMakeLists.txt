cmake_minimum_required(VERSION 3.16)
project(particle_filter)
set_property(GLOBAL PROPERTY pf_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)

# add libraries
add_library(ParticleFilter SHARED src/detail/particle_filter.cpp)
target_include_directories(ParticleFilter PUBLIC include ${EIGEN3_INCLUDE_DIR})
# target_link_libraries(ParticleFilter PUBLIC ginac cln)


# add examples
add_executable(pf_test tests/pf_test.cpp)
target_include_directories(pf_test PUBLIC include)
target_link_libraries(pf_test PUBLIC ParticleFilter)
