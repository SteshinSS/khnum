cmake_minimum_required(VERSION 3.12)
project(CFlex)

set(CMAKE_CXX_STANDARD 17)

add_executable(CFlex src/main.cpp src/utilities/reaction_struct.h src/interface/cli.cpp src/interface/cli.h src/parser/model_parser.cpp src/parser/model_parser.h src/modeller/create_metabolite_list.cpp src/modeller/create_metabolite_list.h src/modeller/create_stoichiometry_matrix.cpp src/modeller/create_stoichiometry_matrix.h src/parser/model_supplementary_parser.cpp src/parser/model_supplementary_parser.h src/utilities/measurements_struct.h src/parser/parser_utilities.h src/math/math_utilites.h src/math/flux_balance_analysis.cpp src/math/flux_balance_analysis.h src/utilities/linear_problem.h src/utilities/EMU.h src/utilities/EMU.cpp src/modeller/create_emu_reactions.cpp src/modeller/create_emu_reactions.h src/utilities/reaction_struct.cpp src/utilities/debug_utilites.cpp src/utilities/debug_utilites.h src/utilities/input_substrate.h src/modeller/create_emu_list.cpp src/modeller/create_emu_list.h src/utilities/MID.h src/math/calculate_input_mid.cpp src/math/calculate_input_mid.h src/modeller/create_emu_networks.cpp src/modeller/create_emu_networks.h)

target_link_libraries(CFlex /usr/local/lib/libglpk.a)