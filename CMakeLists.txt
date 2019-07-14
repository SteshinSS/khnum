cmake_minimum_required(VERSION 3.12)
project(CFlex)

file(GLOB_RECURSE sources      src/*.cpp src/*.h)
file(GLOB_RECURSE test test/*.cpp)
file(GLOB_RECURSE libs lib/*.cpp lib/*.h)

add_executable(CFlex ${sources} ${libs} src/interface/file_system.cpp src/interface/file_system.h src/clusterizer/clusterizer.cpp src/clusterizer/clusterizer.h src/parser/parser_results.h src/parser/open_flux_parser.cpp src/parser/open_flux_parser.h src/modeller/modeller.cpp src/modeller/modeller.h src/solver/solver.cpp src/solver/solver.h src/simulator/Simulator.cpp src/simulator/Simulator.h src/solver/solver_package.h)

target_compile_options(CFlex PUBLIC -std=c++17)

target_include_directories(CFlex PUBLIC src src/parser src/solver src/modeller src/interface src/utilities lib/alglib lib/Eigen)


