cmake_minimum_required(VERSION 3.12)
project(CFlex)

file(GLOB_RECURSE sources      src/*.cpp src/*.h)
file(GLOB_RECURSE test test/*.cpp)
file(GLOB_RECURSE libs lib/*.cpp lib/*.h)

add_executable(CFlex ${sources} ${libs} src/math/mfa_math.h src/math/calculate_flux_bounds.cpp src/math/calculate_flux_bounds.h)

target_compile_options(CFlex PUBLIC -std=c++17)

target_include_directories(CFlex PUBLIC src src/parser src/math src/modeller src/interface src/utilities lib/alglib lib/Eigen lib/glpk)

target_link_libraries(CFlex /usr/local/lib/libglpk.a)

