#pragma once

#include "alglib/ap.h"
#include <random>

#include "alglib/optimization.h"
#include "utilities/problem.h"
#include "simulator/simulator.h"

// See http://www.alglib.net/optimization/levenbergmarquardt.php
// and http://www.alglib.net/translator/man/manual.cpp.html#example_minlm_d_v before reading this code

class Solver {
public:
    static Solver *getSolver(const Problem &problem);

    std::vector<alglib::real_1d_array> getResult();

    void Solve();

private:
    static Solver *instance;

    Solver(const Problem &problem);

    void FillBoundVectors();

    void SetOptimizationParameters();

    void GenerateInitialPoints(std::mt19937 &random_source);

    void PrintStartMessage();

    alglib::real_1d_array RunOptimization();

    double GetSSR(const alglib::real_1d_array &residuals);

    void PrintFinalMessage(const alglib::real_1d_array &free_fluxes);

    static std::vector<Flux> CalculateAllFluxesFromFree(const alglib::real_1d_array &free_fluxes_alglib);

    static void Fillf0Array(alglib::real_1d_array &residuals, const std::vector<EmuAndMid> &simulated_mids);

    static void CalculateResidual(const alglib::real_1d_array &free_fluxes,
                                  alglib::real_1d_array &residuals, void *ptr);

private:
    static int iteration_total_;
    static int iteration_;
    static int nullity_;
    static int reactions_num_;
    static int measurements_count_;

    static std::vector<Reaction> reactions_;
    static std::vector<Emu> measured_isotopes_;
    static Matrix nullspace_;
    static std::vector<EmuNetwork> networks_;
    static std::vector<EmuAndMid> input_mids_;
    static std::vector<Measurement> measurements_;

    static alglib::real_1d_array free_fluxes_;
    static alglib::real_1d_array lower_bounds_;
    static alglib::real_1d_array upper_bounds_;

    static alglib::minlmstate state;
    static alglib::minlmreport report;

    static std::vector<alglib::real_1d_array> allSolutions;

};
