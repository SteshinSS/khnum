#include "interface/cli.h"

#include <iostream>
#include <exception>
#include <vector>
#include <memory>
#include <thread>
#include "alglib/ap.h"
#include <chrono>

#include "simulator/generator.h"
#include "modeller/modeller.h"
#include "parser/open_flux_parser/open_flux_parser.h"
#include "solver/solver.h"
#include "clusterizer/clusterizer.h"
#include "parser/python_parser.h"


namespace khnum {
void RunCli() {
    try {
        //std::unique_ptr<IParser> parser(new ParserOpenFlux("../modelLast"));
        std::unique_ptr<IParser> parser(new ParserMaranas());
        parser->ParseExcludedMetabolites();
        parser->ParseMeasuredIsotopes();
        parser->ParseMeasurements();
        parser->ParseCorrectionMatrices();
        parser->ParseReactions();
        parser->ParseSubstrateInput();

        Modeller modeller(parser->GetResults());

        modeller.CalculateInputSubstrateMids();
        modeller.CreateEmuNetworks();
        modeller.CreateNullspaceMatrix();
        modeller.CalculateFluxBounds();
        modeller.CalculateMeasurementsCount();
        modeller.CheckModelForErrors();

        Problem problem = modeller.GetProblem();
        SimulatorGenerator generator(problem.simulator_parameters_);
        std::vector<alglib::real_1d_array> allSolutions;

        bool use_multithread = false;
        if (use_multithread) {
            const unsigned int num_threads = std::thread::hardware_concurrency();
            std::cout << num_threads << std::endl;

            std::vector<std::vector<alglib::real_1d_array>> one_thread_solutions(num_threads);
            auto one_thread = [&generator](const Problem &problem, std::vector<alglib::real_1d_array> &result) {
                Solver solver(problem, generator);
                solver.Solve();
                result = solver.GetResult();
            };

            std::vector<std::thread> threads;
            for (int i = 0; i < num_threads; ++i) {
                threads.push_back(std::thread(one_thread, std::ref(problem), std::ref(one_thread_solutions[i])));
            }
            std::chrono::time_point<std::chrono::system_clock> start, end;
            start = std::chrono::system_clock::now();
            for (std::thread &t : threads) {
                t.join();
            }
            end = std::chrono::system_clock::now();
            double elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>
                    (end-start).count();

            elapsed_milliseconds /= 1000;
            std::cout << "Average time: " << static_cast<double>(elapsed_milliseconds) / (num_threads * 30)
                << " seconds per iteration" << std::endl;

            for (const auto &vec : one_thread_solutions) {
                for (const alglib::real_1d_array &solution : vec) {
                    allSolutions.push_back(solution);
                }
            }

        } else {
            Solver solver(problem, generator);
            solver.Solve();
            allSolutions = solver.GetResult();
        }

        Clasterizer clusterizer(allSolutions);
        clusterizer.Start();

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
} // namespace khnum