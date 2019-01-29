#include "metabolic_flux_analysis.h"

#include "../utilities/reaction_struct.h"
#include "../utilities/EMU.h"
#include "../utilities/MID.h"
#include "../math/calculate_mids.h"
#include "../Eigen/Eigen"
#include "../math/math_utilites.h"

#include "../alglib/stdafx.h"
#include "../alglib/optimization.h"

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <random>
#include <limits>

std::map<std::string, Flux> EstimateFluxes(ObjectiveParameters *parameters,
                                           const std::map<std::string, FluxVariability> &flux_ranges,
                                           const Matrix &stoichiometry_matrix,
                                           const std::vector<Reaction> &reactions,
                                           const int iteration_total) {
    // find nullspace
    Eigen::FullPivLU<Matrix> lu_decomposition(stoichiometry_matrix);
    Matrix nullspace = lu_decomposition.kernel();
    parameters->nullspace = &nullspace;

    const int measurements_count = GetMeasurementsCount(parameters);

    // find free fluxes
    const int fluxes_total = stoichiometry_matrix.cols();
    const int nullity = fluxes_total - lu_decomposition.rank();
    alglib::real_1d_array free_fluxes;
    alglib::real_1d_array lower_bounds;
    alglib::real_1d_array upper_bounds;

    free_fluxes.setlength(nullity);
    lower_bounds.setlength(nullity);
    upper_bounds.setlength(nullity);

    std::random_device randomizator;
    std::mt19937 random_source(randomizator());

    const double epsx = 0.00000000001;

    FillBoundVectors(lower_bounds, upper_bounds,
                     flux_ranges, reactions, nullity);

    GenerateInitialPoints(free_fluxes, lower_bounds, upper_bounds, reactions, nullity, random_source);



    alglib::ae_int_t maxits = 0;
    alglib::minlmstate state;
    alglib::minlmreport report;

    alglib::minlmcreatev(nullity, measurements_count, free_fluxes, 0.0000000001, state);
    alglib::minlmsetcond(state, epsx, maxits);
    alglib::minlmsetbc(state, lower_bounds, upper_bounds);

    alglib::minlmoptimize(state, CalculateResidual, NULL, parameters);

    alglib::real_1d_array best_free_fluxes;
    alglib::minlmresults(state, best_free_fluxes, report);

    std::map<std::string, Flux> best_all_fluxes = CalculateAllFluxesFromFree(
            free_fluxes, nullspace, *(parameters->reactions));


    std::vector<EMUandMID> best_simulated_mids = CalculateMids(best_all_fluxes,
                                                          *(parameters->networks),
                                                          *(parameters->input_mids),
                                                          *(parameters->measured_isotopes));

    alglib::real_1d_array best_residuals;
    best_residuals.setlength(measurements_count);

    Fillf0Array(best_residuals, best_simulated_mids, *parameters);

    double best_ssr = GetSSR(best_residuals, measurements_count);

    std::cerr << "0 iteration. SSR: " << best_ssr << std::endl;

    for (int iteration = 1; iteration < iteration_total; ++iteration) {
        GenerateInitialPoints(free_fluxes, lower_bounds, upper_bounds, reactions, nullity, random_source);
        alglib::minlmrestartfrom(state, free_fluxes);
        alglib::minlmoptimize(state, CalculateResidual, NULL, parameters);
        alglib::real_1d_array new_free_fluxes;
        alglib::minlmresults(state, new_free_fluxes, report);

        // check which results are better
        std::map<std::string, Flux> new_all_fluxes = CalculateAllFluxesFromFree(
                free_fluxes, nullspace, *(parameters->reactions));


        std::vector<EMUandMID> new_simulated_mids = CalculateMids(new_all_fluxes,
                                                                   *(parameters->networks),
                                                                   *(parameters->input_mids),
                                                                   *(parameters->measured_isotopes));

        alglib::real_1d_array new_residuals;
        new_residuals.setlength(measurements_count);

        Fillf0Array(new_residuals, new_simulated_mids, *parameters);

        double new_ssr = GetSSR(new_residuals, measurements_count);

        std::cerr << iteration << " iteration. SSR: " << new_ssr << std::endl;

        if (new_ssr < best_ssr) {
            best_ssr = new_ssr;
            best_all_fluxes = new_all_fluxes;
        }
    }

    return best_all_fluxes;
}

void CalculateResidual(const alglib::real_1d_array &free_fluxes,
                       alglib::real_1d_array &residuals,
                       void *ptr) {
    ObjectiveParameters &parameters = *(static_cast<ObjectiveParameters *>(ptr));
    Matrix nullspace = *(parameters.nullspace);

    std::map<std::string, Flux> calculated_fluxes = CalculateAllFluxesFromFree(
            free_fluxes, nullspace, *(parameters.reactions));


    std::vector<EMUandMID> simulated_mids = CalculateMids(calculated_fluxes,
                                                          *(parameters.networks),
                                                          *(parameters.input_mids),
                                                          *(parameters.measured_isotopes));

    Fillf0Array(residuals, simulated_mids, parameters);
    return;
}


void FillBoundVectors(alglib::real_1d_array &lower_bounds,
                      alglib::real_1d_array &upper_bounds,
                      const std::map<std::string, FluxVariability> &flux_ranges,
                      const std::vector<Reaction> &reactions,
                      const int nullity) {
    const int reaction_total = reactions.size();

    for (int i = 0; i < nullity; ++i) {
        lower_bounds[i] = (flux_ranges.at(reactions[reaction_total - nullity + i].name)).lower_bound;
        upper_bounds[i] = (flux_ranges.at(reactions[reaction_total - nullity + i].name)).upper_bound;

    }

    return;
}

void GenerateInitialPoints(alglib::real_1d_array &free_fluxes,
                           const alglib::real_1d_array &lower_bounds,
                           const alglib::real_1d_array &upper_bounds,
                           const std::vector<Reaction> &reactions,
                           const int nullity,
                           std::mt19937 &random_source) {

    const int reaction_total = reactions.size();
    std::uniform_real_distribution<> get_random_point(0.0, 1.0);

    for (int i = 0; i < nullity; ++i) {
        free_fluxes[i] = lower_bounds[i] + get_random_point(random_source) * (upper_bounds[i] - lower_bounds[i]);

    }

    return;
}

int GetMeasurementsCount(ObjectiveParameters *parameters) {
    int measurements_count = 0;

    for (const Measurement &measurement : *(parameters->measurements)) {
        measurements_count += measurement.mid.size();
    }
    return measurements_count;
}

std::map<std::string, Flux> CalculateAllFluxesFromFree(const alglib::real_1d_array &free_fluxes,
                                                       const Matrix &nullspace,
                                                       const std::vector<Reaction> &reactions) {
    Eigen::VectorXd free_fluxes_eigen(free_fluxes.length());
    for (int i = 0; i < free_fluxes.length(); ++i) {
        free_fluxes_eigen[i] = free_fluxes[i];
    }

    Matrix all_fluxes_matrix = nullspace * free_fluxes_eigen;

    std::map<std::string, Flux> all_fluxes;
    for (int i = 0; i < all_fluxes_matrix.rows(); ++i) {
        all_fluxes[reactions.at(i).name] = all_fluxes_matrix(i, 0);
    }

    for (int i = all_fluxes_matrix.rows(); i < reactions.size(); ++i) {
        all_fluxes[reactions.at(i).name] = free_fluxes_eigen(i - all_fluxes_matrix.rows());
    }

    return all_fluxes;
}

void Fillf0Array(alglib::real_1d_array &residuals,
                 const std::vector<EMUandMID> &simulated_mids,
                 const ObjectiveParameters &parameters) {
    int total_residuals = 0;
    for (int isotope = 0; isotope < simulated_mids.size(); ++isotope) {
        for (int mass_shift = 0; mass_shift < simulated_mids[isotope].mid.size(); ++mass_shift) {
            residuals[total_residuals] = simulated_mids[isotope].mid[mass_shift];
            residuals[total_residuals] -= (*(parameters.measurements))[isotope].mid[mass_shift];
            residuals[total_residuals] /= 1 + (*(parameters.measurements))[isotope].errors[mass_shift];
            ++total_residuals;
        }
    }
}

double GetSSR(const alglib::real_1d_array &residuals, int measurements_count) {
    double answer = 0.0;
    for (int measurement = 0; measurement < measurements_count; ++measurement) {
        answer += residuals(measurement);
    }
    return answer;
}

