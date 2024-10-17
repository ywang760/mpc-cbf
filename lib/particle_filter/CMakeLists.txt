cmake_minimum_required(VERSION 3.16)
project(particle_filter)
set_property(GLOBAL PROPERTY pf_ADDED "ON")

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

# add libraries
add_library(ParticleFilter SHARED src/detail/particle_filter.cpp)
target_include_directories(ParticleFilter PUBLIC include)
# target_link_libraries(ParticleFilter PUBLIC ginac cln)


# add examples
add_executable(pf_test tests/pf_test.cpp)
target_include_directories(pf_test PUBLIC include)
target_link_libraries(pf_test PUBLIC ParticleFilter)
