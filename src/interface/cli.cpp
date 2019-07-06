#include "cli.h"

#include "file_system.h"

#include "parser.h"
#include "modeller.h"
#include "mfa_math.h"
#include "utilities.h"

#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <string>


void RunCli() {
    try {
        FileSystem model("../modelTca");

        std::vector<Reaction> reactions = ParseReactions(model.getModelFile());
        reactions = SortReactionsByType(reactions);
        std::vector<EMU> measured_emus = ParseMeasuredIsotopes(model.getMeasuredIsotopesFile());
        std::vector<Measurement> measurements = ParseMeasurments(model.getMeasurementsFile(), measured_emus);
        std::vector<InputSubstrate> input_substrates = ParseInputSubstrates(model.getSubstrateInputFile());


        // Creates all elementary reaction in term of EMU
        std::vector<EMUReaction> all_emu_reactions = CreateAllEMUReactions(reactions, measured_emus);

        // Creates all EMUs of input substrate which take part in all_emu_reactions
        std::vector<EMU> input_emu_list = CreateInputEMUList(all_emu_reactions, input_substrates);

        // Calculates MID vector of every input emu of the input_emu_list
        std::vector<EMUandMID> input_substrates_mids = CalculateInputMid(input_substrates, input_emu_list);

        // Create EMU networks. See Antoniewicz 2007
        std::vector<EMUNetwork> emu_networks = CreateEMUNetworks(all_emu_reactions, input_emu_list, measured_emus);





        // Creates list with all reactions metabolites
        std::vector<std::string> full_metabolite_list = CreateFullMetaboliteList(reactions);

        std::vector<std::string> excluded_metabolites = ParseExcludedMetabolites(model.getExcludedMetabolitesFile());
        std::vector<std::string> included_metabolites = CreateIncludedMetaboliteList(
                full_metabolite_list, excluded_metabolites);


        Matrix stoichiometry_matrix = CreateStoichiometryMatrix(reactions, included_metabolites);

        Matrix nullspace = GetNullspace(stoichiometry_matrix);

        reactions = CalculateFluxBounds(reactions);

        // Pack parameters
        ObjectiveParameters parameters;
        parameters.measured_isotopes = &measured_emus;
        parameters.networks = &emu_networks;
        parameters.reactions = &reactions;
        parameters.input_mids = &input_substrates_mids;
        parameters.measurements = &measurements;
        parameters.nullspace = &nullspace;

        std::vector<Flux> answer = EstimateFluxes(parameters, 10);

        reactions = SortReactionByID(reactions);

        for (const Reaction &reaction : reactions) {
            std::cout << reaction.name << " " << answer[reaction.id] << std::endl;
        }


    } catch (std::runtime_error &error) {
        std::cerr << error.what() << std::endl;
    }

    return;
}
