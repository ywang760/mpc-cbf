cmake_minimum_required(VERSION 3.16)
project(cbf)
set_property(GLOBAL PROPERTY cbf_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)

# add local library


# add libraries
add_library(FovCBF SHARED src/cbf.cpp)
target_include_directories(FovCBF PUBLIC include)
target_link_libraries(FovCBF PUBLIC ginac cln)


# add examples
add_executable(cbf_test tests/cbf_test.cpp)
target_include_directories(cbf_test PUBLIC include)
target_link_libraries(cbf_test PUBLIC ginac cln)

add_executable(test_class tests/test_class.cpp)
target_include_directories(test_class PUBLIC include)
target_link_libraries(test_class PUBLIC FovCBF ginac cln)


