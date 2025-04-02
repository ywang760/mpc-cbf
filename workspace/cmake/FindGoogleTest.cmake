# FindGoogleTest.cmake
# Find the GoogleTest library and configure it for use

# Define search locations - using relative paths from this cmake file
set(GOOGLETEST_SEARCH_PATHS
    ${CMAKE_CURRENT_LIST_DIR}/../../third_party/googletest  # From cmake dir to project root and then to third_party
    ${CMAKE_SOURCE_DIR}/third_party/googletest  # From top-level project
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/googletest  # Alternative location
    ${CMAKE_SOURCE_DIR}/../third_party/googletest  # One level up
)

# Try to find GoogleTest
set(GOOGLETEST_FOUND FALSE)

foreach(SEARCH_PATH IN LISTS GOOGLETEST_SEARCH_PATHS)
    if(EXISTS ${SEARCH_PATH} AND NOT GOOGLETEST_FOUND)
        set(GOOGLETEST_ROOT ${SEARCH_PATH})
        set(GOOGLETEST_FOUND TRUE)
        break()
    endif()
endforeach()

# Configure GoogleTest if found
if(GOOGLETEST_FOUND)
    # Only add GoogleTest if it hasn't been added yet
    if(NOT TARGET gtest AND NOT TARGET gmock)
        message(STATUS "Configuring GoogleTest from ${GOOGLETEST_ROOT}")
        # Make sure both gtest and gmock are built
        set(BUILD_GMOCK ON CACHE BOOL "Build GMock" FORCE)
        add_subdirectory(${GOOGLETEST_ROOT} ${CMAKE_BINARY_DIR}/googletest EXCLUDE_FROM_ALL)
    endif()

    # Set variables to indicate that GoogleTest was found
    set(GOOGLETEST_CONFIGURED TRUE)
    set(GMOCK_INCLUDE_DIRS ${GOOGLETEST_ROOT}/googlemock/include)
    set(GTEST_INCLUDE_DIRS ${GOOGLETEST_ROOT}/googletest/include)

    include(CTest)
    include(GoogleTest)
else()
    message(WARNING "GoogleTest not found in any of the search paths")
    set(GOOGLETEST_CONFIGURED FALSE)
endif()
