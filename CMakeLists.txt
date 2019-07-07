cmake_minimum_required(VERSION 3.12)
project(CFlex)

file(GLOB_RECURSE sources      src/*.cpp src/*.h)
file(GLOB_RECURSE test test/*.cpp)
file(GLOB_RECURSE libs lib/*.cpp lib/*.h)

add_executable(CFlex ${sources} ${libs} src/math/mfa_math.h src/math/calculate_flux_bounds.cpp src/math/calculate_flux_bounds.h src/interface/file_system.cpp src/interface/file_system.h src/math/clusterizer.cpp src/math/clusterizer.h src/parser/parser_results.h src/parser/parser_open_flux.cpp src/parser/parser_open_flux.h)

target_compile_options(CFlex PUBLIC -std=c++17)

target_include_directories(CFlex PUBLIC src src/parser src/math src/modeller src/interface src/utilities lib/alglib lib/Eigen)


