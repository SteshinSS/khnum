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
#include "../math/metabolic_flux_analysis.h"
#include "../modeller/sort_reactions.h"
#include "../utilities/objective_parameters.h"

#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <string>


void RunCli() {
    try {


        // ../parser/model_parser.h
        std::vector<Reaction> reactions = ParseReactions("../model/model.csv");

        // ../modeller/sort_reactions.h
        reactions = SortReactions(reactions);

        // ../parser/model_supplementary_parser.h
        std::vector<EMU> measured_emus = ParseMeasuredIsotopes("../model/measured_isotopes.txt");

        // ../parser/model_supplementary_parser.h
        std::vector<Measurement> measurements = ParseMeasurments("../model/measurements.csv", measured_emus);

        // ../modeller/create_emu_reactions.h
        // Creates all elementary reaction in term of EMU
        std::vector<EMUReaction> all_emu_reactions = CreateAllEMUReactions(reactions, measured_emus);

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
        std::vector<EMUNetwork> emu_networks = CreateEMUNetworks(all_emu_reactions, input_emu_list, measured_emus);




        // ../modeller/create_metabolite_list.h
        std::vector<std::string> full_metabolite_list = CreateFullMetaboliteList(reactions);

        // ../parser/model_supplementary_parser.h
        std::vector<std::string> excluded_metabolites = ParseExcludedMetabolites("../model/excluded_metabolites.txt");
        std::vector<std::string> included_metabolites = CreateIncludedMetaboliteList(
                full_metabolite_list, excluded_metabolites);

        // ../modeller/create_stoichiometry_matrix.h
        Matrix stoichiometry_matrix = CreateStoichiometryMatrix(reactions, included_metabolites);

        std::cerr << stoichiometry_matrix << std::endl;

        // ../math/flux_balance_analysis.h
        // Run preliminary FBA for calculate initial fluxes
        // std::map<std::string, Flux> initial_fluxes = EstablishInitialFluxes(
        //        stoichiometry_matrix, reactions, included_metabolites);

        // std::vector<FluxVariability> flux_ranges = EstablishAllFluxRanges(
        //        stoichiometry_matrix, reactions, included_metabolites);

        for (Reaction &reaction : reactions) {
            if (std::isnan(reaction.basis)) {
                reaction.computed_upper_bound = reaction.setted_upper_bound;
                reaction.computed_lower_bound = reaction.setted_lower_bound;
            } else {
                if (std::isnan(reaction.deviation)) {
                    reaction.computed_upper_bound = reaction.basis;
                    reaction.computed_lower_bound = reaction.basis;
                } else {
                    reaction.computed_upper_bound = reaction.basis + reaction.deviation;
                    reaction.computed_lower_bound = reaction.basis - reaction.deviation;
                }
            }
        }

        // Pack parameters
        ObjectiveParameters parameters;
        parameters.measured_isotopes = &measured_emus;
        parameters.networks = &emu_networks;
        parameters.reactions = &reactions;
        parameters.input_mids = &input_substrates_mids;
        parameters.measurements = &measurements;
        parameters.nullspace = nullptr;

        std::vector<Flux> answer = EstimateFluxes(&parameters,
                                                  stoichiometry_matrix,
                                                  reactions,
                                                  10);

        for (const Reaction &reaction : reactions) {
            std::cerr << reaction.name << " " << answer[reaction.id] << std::endl;
        }


    } catch (std::runtime_error &error) {
        std::cerr << error.what() << std::endl;
    }

    return;
}