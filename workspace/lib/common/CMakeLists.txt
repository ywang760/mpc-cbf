cmake_minimum_required(VERSION 3.21)

project(common)
set_property(GLOBAL PROPERTY common_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../cmake ${CMAKE_MODULE_PATH})

if (NOT TARGET spdlog)
    add_subdirectory(../../../third_party/spdlog
        ${CMAKE_CURRENT_BINARY_DIR}/spdlog)
endif()

# find packages
find_package(Eigen3 REQUIRED)

# add libraries
add_library(common_Logging INTERFACE)
target_include_directories(common_Logging INTERFACE
    ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(common_Logging INTERFACE
    Eigen3::Eigen
    ginac
    spdlog::spdlog_header_only
)