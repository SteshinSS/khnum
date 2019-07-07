#include "cli.h"

#include "file_system.h"

#include "parser_open_flux.h"
#include "modeller.h"
#include "mfa_math.h"
#include "utilities.h"

#include "clusterizer.h"

#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <string>
#include <memory>


void RunCli() {
    try {
        FileSystem model("../modelTca");

        /*
        std::vector<Reaction> reactions = ParseReactions(model.getModelFile());
        reactions = SortReactionsByType(reactions);
        std::vector<Emu> measured_emus = ParseMeasuredIsotopes(model.getMeasuredIsotopesFile());
        std::vector<Measurement> measurements = ParseMeasurments(model.getMeasurementsFile(), measured_emus);
        std::vector<InputSubstrate> input_substrates = ParseInputSubstrates(model.getSubstrateInputFile()); */

        std::unique_ptr<Parser> parser(new ParserOpenFlux("../modelTca"));
        parser->ReadExcludedMetabolites();
        parser->ReadMeasuredIsotopes();
        parser->ReadMeasurements();
        parser->ReadReactions();
        parser->ReadSubstrateInput();


        // Creates all elementary reaction in term of Emu
        std::vector<EMUReaction> all_emu_reactions = CreateAllEMUReactions(*parser->GetResults().reactions,
                                                                           *parser->GetResults().measuredEmu);

        // Creates all EMUs of input substrate which take part in all_emu_reactions
        std::vector<Emu> input_emu_list = CreateInputEMUList(all_emu_reactions, *parser->GetResults().input_substrate);

        // Calculates MID vector of every input emu of the input_emu_list
        std::vector<EmuAndMid> input_substrates_mids = CalculateInputMid(*parser->GetResults().input_substrate,
                                                                         input_emu_list);

        // Create Emu networks. See Antoniewicz 2007
        std::vector<EMUNetwork> emu_networks = CreateEMUNetworks(all_emu_reactions, input_emu_list,
                                                                 *parser->GetResults().measuredEmu);





        // Creates list with all reactions metabolites
        std::vector<std::string> full_metabolite_list = CreateFullMetaboliteList(*parser->GetResults().reactions);

        //std::vector<std::string> excluded_metabolites = ParseExcludedMetabolites(model.getExcludedMetabolitesFile());
        std::vector<std::string> included_metabolites = CreateIncludedMetaboliteList(
                full_metabolite_list, *parser->GetResults().excluded_metabolites);


        Matrix stoichiometry_matrix = CreateStoichiometryMatrix(*parser->GetResults().reactions, included_metabolites);

        Matrix nullspace = GetNullspace(stoichiometry_matrix);

        auto reactions = CalculateFluxBounds(*parser->GetResults().reactions);

        auto pack = parser->GetResults();
        // Pack parameters
        ObjectiveParameters parameters;
        parameters.measured_isotopes = &(*pack.measuredEmu);
        parameters.networks = &emu_networks;
        parameters.reactions = &reactions;
        parameters.input_mids = &input_substrates_mids;
        parameters.measurements = &(*pack.measurements);
        parameters.nullspace = &nullspace;

        std::vector<alglib::real_1d_array> allSolutions = EstimateFluxes(parameters, 10);

        Clasterizer a(allSolutions);
/*
        reactions = SortReactionByID(reactions);

        for (const Reaction &reaction : reactions) {
            std::cout << reaction.name << " " << answer[reaction.id] << std::endl;
        }
*/

    } catch (std::runtime_error &error) {
        std::cerr << error.what() << std::endl;
    }

    return;
}
