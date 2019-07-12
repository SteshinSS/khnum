#include "cli.h"

#include "file_system.h"

#include "modeller/modeller.h"
#include "parser/open_flux_parser.h"
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
        std::unique_ptr<Parser> parser(new ParserOpenFlux("../modelTca"));
        parser->ReadExcludedMetabolites();
        parser->ReadMeasuredIsotopes();
        parser->ReadMeasurements();
        parser->ReadReactions();
        parser->ReadSubstrateInput();

        Modeller modeller(parser->GetResults());

        modeller.CalculateInputSubstrateMids();
        modeller.CreateEmuNetworks();
        modeller.CreateNullspaceMatrix();
        modeller.CalculateFluxBounds();

        Problem problem = modeller.GetProblem();

        std::vector<alglib::real_1d_array> allSolutions = EstimateFluxes(problem, 10);

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
