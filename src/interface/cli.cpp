#include "cli.h"
#include "../parser/model_parser.h"
#include "../parser/model_supplementary_parser.h"
#include "../modeller/create_metabolite_list.h"
#include "../utilities/EMU.h"
#include "../modeller/create_emu_reactions.h"
#include "../modeller/create_stoichiometry_matrix.h"
#include "../modeller/create_emu_networks.h"
#include "../math/math_utilites.h"
#include "../math/flux_balance_analysis.h"
#include "../utilities/debug_utilites.h"
#include "../utilities/input_substrate.h"
#include "../modeller/create_emu_list.h"
#include "../math/calculate_input_mid.h"
#include "../utilities/MID.h"
#include "../math/calculate_mids.h"

#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <string>

void RunCli() {
    try {

        // ../parser/model_parser.h
        std::vector<Reaction> reactions = ParseReactions("../model/model.csv");

        // ../parser/model_supplementary_parser.h
        std::vector<EMU> measured_isotopes = ParseMeasuredIsotopes("../model/measured_isotopes.txt");

        // ../modeller/create_emu_reactions.h
        // Creates all elementary reaction in term of EMU
        std::vector<EMUReaction> all_emu_reactions = CreateAllEMUReactions(reactions, measured_isotopes);

        // ../parser/model_supllemetary_parser.h
        std::vector<InputSubstrate> input_substrates = ParseInputSubstrates("../model/substrate_input.csv");

        // ../modeller/create_emu_list.h
        // Creates all EMUs of input substrate which take part in all_emu_reactions
        std::vector<EMU> input_emu_list = CreateInputEMUList(all_emu_reactions, input_substrates);

        // ../math/calculate_input_mid.h
        // Calculate MID vector of every input emu from the input_emu_list
        std::vector<EMUandMID> input_substrates_mids = CalculateInputMid(input_substrates, input_emu_list);

        // ../modeller/create_emu_networks.h
        // Create EMU networks. See Antoniewicz 2007
        std::vector<EMUNetwork> emu_networks = CreateEMUNetworks(all_emu_reactions, input_emu_list, measured_isotopes);




        // ../modeller/create_metabolite_list.h
        std::vector<std::string> full_metabolite_list = CreateFullMetaboliteList(reactions);

        // ../parser/model_supplementary_parser.h
        std::vector<std::string> excluded_metabolites = ParseExcludedMetabolites("../model/excluded_metabolites.txt");
        std::vector<std::string> included_metabolites = CreateIncludedMetaboliteList(
            full_metabolite_list, excluded_metabolites);

        // ../modeller/create_stoichiometry_matrix.h
        Matrix stoichiometry_matrix = CreateStoichiometryMatrix(reactions, included_metabolites);

        // ../math/flux_balance_analysis.h
        // Run preliminary FBA for calculate initial fluxes
        std::map<std::string, Flux> initial_fluxes = EstablishInitialFluxes(
            stoichiometry_matrix, reactions, included_metabolites);

        std::map<std::string, FluxVariability>  flux_ranges = EstablishAllFluxRanges(
            stoichiometry_matrix, reactions, included_metabolites);

        SolveOneNetwork(initial_fluxes, emu_networks[0], input_substrates_mids, 1);


    } catch (std::runtime_error &error) {
        std::cerr << error.what() << std::endl;
    }

    return;
}