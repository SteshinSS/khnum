cmake_minimum_required(VERSION 3.0)
project(khnum)
set(CMAKE_CXX_STANDARD 17)


if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "FATAL: In-source builds are not allowed. You should create a separate directory for build files.")
endif()

file(GLOB_RECURSE khnum_sources src/*.cpp src/*.h)
add_executable(khnum ${khnum_sources})

if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Og -ggdb")
ELSE()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
ENDIF()

target_include_directories(khnum PUBLIC src)
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS "-Wall -Wpedantic -Wextra")


# Build libraries

add_subdirectory(lib/alglib)
add_subdirectory(lib/eigen)
target_link_libraries(khnum eigen alglib)




