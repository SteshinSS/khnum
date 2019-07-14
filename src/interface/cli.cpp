#include "interface/cli.h"

#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <string>
#include <memory>
#include "alglib/ap.h"

#include "modeller/modeller.h"
#include "parser/open_flux_parser.h"
#include "solver/solver.h"
#include "clusterizer/clusterizer.h"


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
        modeller.CalculateMeasurementsCount();

        Problem problem = modeller.GetProblem();

        Solver *solver = Solver::getSolver(problem);
        solver->Solve();
        std::vector<alglib::real_1d_array> allSolutions = solver->getResult();

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
