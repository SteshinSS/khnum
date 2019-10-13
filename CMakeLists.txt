cmake_minimum_required(VERSION 3.0)
project(khnum)
set(CMAKE_CXX_STANDARD 17)


if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "FATAL: In-source builds are not allowed. You should create a separate directory for build files.")
endif()


file(GLOB_RECURSE KHNUM_SOURCES
     src/clusterizer/*
     src/interface/*
     src/modeller/*
     src/parser/*
     src/simulator/*
     src/solver/*
     src/utilities/*)

add_executable(khnum ${KHNUM_SOURCES} src/main.cpp)

if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Og -ggdb")
ELSE()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
ENDIF()

target_include_directories(khnum PUBLIC src)
set_source_files_properties(${KHNUM_SOURCES} PROPERTIES COMPILE_FLAGS "-Wall -Wpedantic -Wextra")


# Build libraries

add_subdirectory(lib/alglib)
add_subdirectory(lib/eigen)
target_link_libraries(khnum eigen alglib)

option(BUILD_TESTS "Determines whether to build tests." ON)
if (BUILD_TESTS)
    enable_testing()
    file(GLOB_RECURSE TESTS_SOURCES tests/*.cpp tests/*.h)

    add_executable(run_tests ${TESTS_SOURCES} ${KHNUM_SOURCES})
    add_subdirectory(lib/catch)
    target_link_libraries(run_tests eigen alglib catch)
    target_include_directories(run_tests PUBLIC src)
    add_test(NAME run_tests.cpp COMMAND MainTests)
endif()