# Test executables for cbf library components

if (NOT TARGET spdlog)
    add_subdirectory(../../../../third_party/spdlog
        ${CMAKE_CURRENT_BINARY_DIR}/spdlog)
endif()

# Configure tests for cbf library
# TODO: add more test sources as needed
set(TEST_SOURCES
    TestInitConnectivity.cpp
    TestInitSafetyCBF.cpp
)

add_executable(cbf_tests ${TEST_SOURCES})

target_compile_definitions(cbf_tests PRIVATE
    SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_DEBUG
)

target_link_libraries(cbf_tests
    PUBLIC
    cbf_detail_ConnectivityCBF
    gtest
    gtest_main
    gmock
    gmock_main
    spdlog::spdlog_header_only
)

# Add include directories
target_include_directories(cbf_tests PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${GTEST_INCLUDE_DIRS}
    ${GMOCK_INCLUDE_DIRS}
)

gtest_discover_tests(cbf_tests)
