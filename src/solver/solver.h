#ifndef CFLEX_SOLVER_H
#define CFLEX_SOLVER_H

#include <ap.h>
#include <random>

#include "optimization.h"
#include "utilities/problem.h"

class Solver {
public:
    static Solver* getInstance(Problem &problem);
    void Solve();
    std::vector<alglib::real_1d_array> getResult();

private:
    static Solver* instance;
    Solver(Problem &problem);

    void FillBoundVectors();
    void SetOptimizationParameters();
    void GenerateInitialPoints();
    void PrintStartMessage();
    alglib::real_1d_array RunOptimization();

    double GetSSR(const alglib::real_1d_array &residuals);


    void PrintFinalMessage(const alglib::real_1d_array &free_fluxes,
                           const double ssr);


    static std::vector<Flux> CalculateAllFluxesFromFree(const alglib::real_1d_array &free_fluxes_alglib);


    static Eigen::VectorXd GetEigenVectorFromAlgLibVector(const alglib::real_1d_array &alglib_vector);


    static void Fillf0Array(alglib::real_1d_array &residuals, const std::vector<EmuAndMid> &simulated_mids);


    static void CalculateResidual(const alglib::real_1d_array &free_fluxes,
                                  alglib::real_1d_array &residuals, void* ptr);


private:
    static int iteration_total_;
    static int iteration_;

    static std::vector<Reaction> reactions_;
    static std::vector<Emu> measured_isotopes_;
    static Matrix nullspace_;
    static std::vector<EMUNetwork> networks_;
    static std::vector<EmuAndMid> input_mids_;
    static std::vector<Measurement> measurements_;
    static int measurements_count_;

    static int nullity_;

    static alglib::real_1d_array free_fluxes_;
    static alglib::real_1d_array lower_bounds_;
    static alglib::real_1d_array upper_bounds_;

    static alglib::minlmstate state;
    static alglib::minlmreport report;

    static std::vector<alglib::real_1d_array> allSolutions;
};

#endif //CFLEX_SOLVER_H
