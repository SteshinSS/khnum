cmake_minimum_required(VERSION 3.12)
project(CFlex)

file(GLOB_RECURSE sources      src/*.cpp src/*.h)
file(GLOB_RECURSE test test/*.cpp)
file(GLOB_RECURSE libs lib/*.cpp lib/*.h)

add_executable(CFlex ${sources} ${libs} src/clusterizer/clusterizer.cpp src/clusterizer/clusterizer.h src/parser/parser_results.h src/parser/open_flux_parser.cpp src/parser/open_flux_parser.h src/modeller/modeller.cpp src/modeller/modeller.h src/solver/solver.cpp src/solver/solver.h src/simulator/simulator.cpp src/simulator/simulator.h src/solver/solver_package.h)

target_compile_options(CFlex PUBLIC -std=c++17)

target_include_directories(CFlex PUBLIC src lib)


