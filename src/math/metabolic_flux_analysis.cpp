#include "metabolic_flux_analysis.h"

#include "reaction_struct.h"
#include "EMU.h"
#include "MID.h"
#include "calculate_mids.h"
#include "Eigen"
#include "math_utilites.h"

#include "stdafx.h"
#include "optimization.h"

#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <tuple>
#include <limits>

std::vector<Flux> EstimateFluxes(ObjectiveParameters &parameters,
                                 const int iteration_total) {

    const int measurements_count = GetMeasurementsCount(*(parameters.measurements));
    const int nullity = (*(parameters.nullspace)).cols();
    const std::vector<Reaction> &reactions = *(parameters.reactions);

    alglib::real_1d_array free_fluxes;
    free_fluxes.setlength(nullity);

    alglib::real_1d_array lower_bounds;
    lower_bounds.setlength(nullity);

    alglib::real_1d_array upper_bounds;
    upper_bounds.setlength(nullity);

    std::random_device randomizator;
    std::mt19937 random_source(randomizator());


    FillBoundVectors(lower_bounds, upper_bounds, reactions, nullity);

    GenerateInitialPoints(free_fluxes, lower_bounds, upper_bounds, reactions, nullity, 1, random_source);

    alglib::minlmstate state;
    alglib::minlmreport report;

    SetOptimizationParameters(free_fluxes, lower_bounds, upper_bounds, nullity, measurements_count, state,
                              report);

    auto[best_ssr, best_fluxes] = RunOptimization(measurements_count,
                                                  &parameters, state, report);

    for (int iteration = 1; iteration < iteration_total; ++iteration) {
        GenerateInitialPoints(free_fluxes, lower_bounds, upper_bounds, reactions, nullity, iteration + 1, random_source);
        alglib::minlmrestartfrom(state, free_fluxes);

        auto[new_ssr, new_fluxes] = RunOptimization(measurements_count,
                                                    &parameters, state, report);

        if (new_ssr < best_ssr) {
            best_ssr = new_ssr;
            best_fluxes = new_fluxes;
        }
    }

    return best_fluxes;
}


int GetMeasurementsCount(const std::vector<Measurement> &measurements) {
    int measurements_count = 0;

    for (const Measurement &measurement : measurements) {
        measurements_count += measurement.mid.size();
    }
    return measurements_count;
}


void FillBoundVectors(alglib::real_1d_array &lower_bounds,
                      alglib::real_1d_array &upper_bounds,
                      const std::vector<Reaction> &reactions,
                      const int nullity) {
    const int reaction_total = reactions.size();

    for (int i = 0; i < nullity; ++i) {
        lower_bounds[i] = reactions[reaction_total - nullity + i].computed_lower_bound;
        upper_bounds[i] = reactions[reaction_total - nullity + i].computed_upper_bound;
    }

}


void GenerateInitialPoints(alglib::real_1d_array &free_fluxes,
                           const alglib::real_1d_array &lower_bounds,
                           const alglib::real_1d_array &upper_bounds,
                           const std::vector<Reaction> &reactions,
                           const int nullity,
                           const int iteration,
                           std::mt19937 &random_source) {

    std::uniform_real_distribution<> get_random_point(0.0, 1.0);

    for (int i = 0; i < nullity; ++i) {
        free_fluxes[i] = lower_bounds[i] + get_random_point(random_source) * (upper_bounds[i] - lower_bounds[i]);
    }

    PrintStartMessage(free_fluxes, reactions, nullity, iteration);
}

void PrintStartMessage(const alglib::real_1d_array &free_fluxes,
                       const std::vector<Reaction> &reactions,
                       const int nullity,
                       const int iteration) {
    std::cout << "Start " << iteration << " iteration from: " << std::endl;
    for (int i = 0; i < nullity; ++i) {
        std::cout << reactions[reactions.size() - nullity + i].name << " = " << free_fluxes[i] << std::endl;
    }
    std::cout << std::endl;
}

void SetOptimizationParameters(alglib::real_1d_array &free_fluxes,
                               alglib::real_1d_array &lower_bounds,
                               alglib::real_1d_array &upper_bounds,
                               const int nullity,
                               const int measurements_count,
                               alglib::minlmstate &state,
                               alglib::minlmreport &report) {
    alglib::ae_int_t maxits = 0;
    const double epsx = 0.00000000001;

    alglib::minlmcreatev(nullity, measurements_count, free_fluxes, 0.0001, state);
    alglib::minlmsetcond(state, epsx, maxits);
    alglib::minlmsetbc(state, lower_bounds, upper_bounds);
}

std::tuple<double, std::vector<Flux>> RunOptimization(int measurements_count,
                                                      ObjectiveParameters *parameters,
                                                      alglib::minlmstate &state,
                                                      alglib::minlmreport &report) {

    alglib::minlmoptimize(state, CalculateResidual, NULL, parameters);

    alglib::real_1d_array final_free_fluxes;
    alglib::minlmresults(state, final_free_fluxes, report);

    std::vector<Flux> final_all_fluxes = CalculateAllFluxesFromFree(
            final_free_fluxes, *(parameters->nullspace), *(parameters->reactions));


    std::vector<EMUandMID> simulated_mids = CalculateMids(final_all_fluxes,
                                                          *(parameters->networks),
                                                          *(parameters->input_mids),
                                                          *(parameters->measured_isotopes));

    alglib::real_1d_array residuals;
    residuals.setlength(measurements_count);

    Fillf0Array(residuals, simulated_mids, *(parameters->measurements));

    double best_ssr = GetSSR(residuals, measurements_count);

    PrintFinalMessage(final_free_fluxes, *(parameters->reactions), best_ssr, report);

    return std::tie(best_ssr, final_all_fluxes);
}


void CalculateResidual(const alglib::real_1d_array &free_fluxes,
                       alglib::real_1d_array &residuals,
                       void *ptr) {
    ObjectiveParameters &parameters = *(static_cast<ObjectiveParameters *>(ptr));

    std::vector<Flux> calculated_fluxes = CalculateAllFluxesFromFree(
            free_fluxes, *(parameters.nullspace), *(parameters.reactions));

    std::vector<EMUandMID> simulated_mids = CalculateMids(calculated_fluxes,
                                                          *(parameters.networks),
                                                          *(parameters.input_mids),
                                                          *(parameters.measured_isotopes));

    Fillf0Array(residuals, simulated_mids, *(parameters.measurements));
}


std::vector<Flux> CalculateAllFluxesFromFree(const alglib::real_1d_array &free_fluxes_alglib,
                                             const Matrix &nullspace,
                                             const std::vector<Reaction> &reactions) {
    Eigen::VectorXd free_fluxes_eigen = GetEigenVectorFromAlgLibVector(free_fluxes_alglib);
    Matrix all_fluxes_matrix = nullspace * free_fluxes_eigen;
    std::vector<Flux> all_fluxes(reactions.size());

    // non metabolite balance reactions
    const int real_reactions_total = all_fluxes_matrix.rows();

    // metabolite balance
    const int fake_reactions_total = reactions.size() - all_fluxes_matrix.rows();

    for (int i = 0; i < real_reactions_total; ++i) {
        all_fluxes[reactions.at(reactions.size() - real_reactions_total + i).id] = all_fluxes_matrix(i, 0);
    }

    // fill const fake fluxes
    for (int i = 0; i < fake_reactions_total; ++i) {
        all_fluxes[reactions.at(i).id] = 1;
    }

    return all_fluxes;
}

Eigen::VectorXd GetEigenVectorFromAlgLibVector(const alglib::real_1d_array &alglib_vector) {
    Eigen::VectorXd eigen_vector(alglib_vector.length());
    for (int i = 0; i < alglib_vector.length(); ++i) {
        eigen_vector(i) = alglib_vector[i];
    }

    return eigen_vector;
}

void Fillf0Array(alglib::real_1d_array &residuals,
                 const std::vector<EMUandMID> &simulated_mids,
                 const std::vector<Measurement> &measurements) {
    int total_residuals = 0;
    for (int isotope = 0; isotope < simulated_mids.size(); ++isotope) {
        for (int mass_shift = 0; mass_shift < simulated_mids[isotope].mid.size(); ++mass_shift) {
            residuals[total_residuals] = simulated_mids[isotope].mid[mass_shift];
            residuals[total_residuals] -= (measurements[isotope].mid[mass_shift]);
            residuals[total_residuals] /= 1 + measurements[isotope].errors[mass_shift];
            ++total_residuals;
        }
    }
}

double GetSSR(const alglib::real_1d_array &residuals, const int measurements_count) {
    double answer = 0.0;
    for (int measurement = 0; measurement < measurements_count; ++measurement) {
        answer += residuals(measurement) * residuals(measurement);
    }
    return answer;
}

void PrintFinalMessage(const alglib::real_1d_array &free_fluxes,
                       const std::vector<Reaction> &reactions,
                       const double ssr,
                       const alglib::minlmreport &report) {
    std::cout << "Finish at: " << std::endl;

    for (int i = 0; i < free_fluxes.length(); ++i) {
        std::cout << reactions[reactions.size() - free_fluxes.length() + i].name <<
                  " = " << free_fluxes[i] << std::endl;
    }
    std::cout << "with SSR: " << ssr << " in " << report.iterationscount << " steps." << std::endl << std::endl;
}