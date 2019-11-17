cmake_minimum_required(VERSION 3.0)
project(khnum)
set(CMAKE_CXX_STANDARD 17)

option(BUILD_TESTS "Determines whether to build tests." OFF)



if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
ELSE()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
ENDIF()

if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "FATAL: In-source builds are not allowed. You should create a separate directory for build files.")
endif()

# Build libraries

add_subdirectory(lib/alglib)
add_subdirectory(lib/eigen)
add_subdirectory(lib/glpk)

# Create khnum library

file(GLOB_RECURSE KHNUM_SOURCES
     src/clusterizer/*
     src/interface/*
     src/modeller/*
     src/parser/*
     src/simulator/*
     src/solver/*
     src/utilities/*)

set_source_files_properties(${KHNUM_SOURCES} PROPERTIES COMPILE_FLAGS "-Wall -Wpedantic -Wextra")

add_library(khnum_lib ${KHNUM_SOURCES})
target_link_libraries(khnum_lib eigen alglib glpk)
target_include_directories(khnum_lib PUBLIC include)

# Create khnum binary

add_executable(khnum src/main.cpp)
target_link_libraries(khnum khnum_lib)
target_include_directories(khnum PUBLIC include)

# Create tests


if (BUILD_TESTS)
    enable_testing()
    file(GLOB_RECURSE TESTS_SOURCES tests/*.cpp tests/*.h)
    set_source_files_properties(${TESTS_SOURCES} PROPERTIES COMPILE_FLAGS "-Og -ggdb")

    add_executable(run_tests ${TESTS_SOURCES})
    target_include_directories(run_tests PUBLIC include)
    add_subdirectory(lib/catch)
    target_link_libraries(run_tests khnum_lib eigen alglib glpk catch)
    add_test(NAME open_flux_parser_test.cpp COMMAND OpenFluxParserTest)
    add_test(NAME open_flux_parse_reactions_test.cpp COMMAND OpenFluxParseReactionsTest)
endif()