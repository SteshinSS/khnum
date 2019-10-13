cmake_minimum_required(VERSION 3.0)
project(khnum_tests)
set(CMAKE_CXX_STANDARD 17)


if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "FATAL: In-source builds are not allowed. You should create a separate directory for build files.")
endif()

file(GLOB_RECURSE test_sources ${TEST_SOURCE_DIR}/src/*.cpp ${TEST_SOURCE_DIR}/src/*.h *.cpp)
add_executable(khnum_tests ${test_sources})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Og -ggdb")

target_include_directories(khnum_tests PUBLIC ${CTEST_SOURCE_DIRECTORY}/src)


# Build libraries

target_link_libraries(khnum_tests eigen alglib catch)