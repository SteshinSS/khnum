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
#include <iostream>
#include <random>
#include <limits>

std::vector<Flux> EstimateFluxes(ObjectiveParameters *parameters,
                                 const Matrix &stoichiometry_matrix,
                                 const std::vector<Reaction> &reactions,
                                 const int iteration_total) {
    Matrix nullspace = GetRREF(stoichiometry_matrix);

    std::cerr << nullspace << std::endl;
    parameters->nullspace = &nullspace;
    const int measurements_count = GetMeasurementsCount(parameters);

    const int nullity = nullspace.cols();
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
                     reactions, nullity);

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

    std::vector<Flux> best_all_fluxes = CalculateAllFluxesFromFree(
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
        std::vector<Flux> new_all_fluxes = CalculateAllFluxesFromFree(
                free_fluxes, nullspace, *(parameters->reactions));


        std::vector<EMUandMID> new_simulated_mids = CalculateMids(new_all_fluxes,
                                                                  *(parameters->networks),
                                                                  *(parameters->input_mids),
                                                                  *(parameters->measured_isotopes));

        alglib::real_1d_array new_residuals;
        new_residuals.setlength(measurements_count);

        Fillf0Array(new_residuals, new_simulated_mids, *parameters);

        double new_ssr = GetSSR(new_residuals, measurements_count);

        std::cerr << iteration << " iteration. SSR: " << new_ssr << " steps: " << report.iterationscount << std::endl;

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

    std::vector<Flux> calculated_fluxes = CalculateAllFluxesFromFree(
            free_fluxes, nullspace, *(parameters.reactions));

    std::vector<EMUandMID> simulated_mids = CalculateMids(calculated_fluxes,
                                                          *(parameters.networks),
                                                          *(parameters.input_mids),
                                                          *(parameters.measured_isotopes));

    Fillf0Array(residuals, simulated_mids, parameters);
}


void FillBoundVectors(alglib::real_1d_array &lower_bounds,
                      alglib::real_1d_array &upper_bounds,
                      const std::vector<Reaction> &reactions,
                      const int nullity) {
    const int reaction_total = reactions.size();

    for (int i = 0; i < nullity; ++i) {
        lower_bounds[i] = reactions[reaction_total - nullity + i].computed_lower_bound;
        upper_bounds[i] = reactions[reaction_total - nullity + i].computed_upper_bound;
        std::cerr << reactions[reaction_total - nullity + i].name << " " << lower_bounds[i] << " " << upper_bounds[i]
                  << std::endl;
    }

}

void GenerateInitialPoints(alglib::real_1d_array &free_fluxes,
                           const alglib::real_1d_array &lower_bounds,
                           const alglib::real_1d_array &upper_bounds,
                           const std::vector<Reaction> &reactions,
                           const int nullity,
                           std::mt19937 &random_source) {

    std::uniform_real_distribution<> get_random_point(0.0, 1.0);

    for (int i = 0; i < nullity; ++i) {
        free_fluxes[i] = lower_bounds[i] + get_random_point(random_source) * (upper_bounds[i] - lower_bounds[i]);
    }


}

int GetMeasurementsCount(ObjectiveParameters *parameters) {
    int measurements_count = 0;

    for (const Measurement &measurement : *(parameters->measurements)) {
        measurements_count += measurement.mid.size();
    }
    return measurements_count;
}

std::vector<Flux> CalculateAllFluxesFromFree(const alglib::real_1d_array &free_fluxes,
                                             const Matrix &nullspace,
                                             const std::vector<Reaction> &reactions) {
    Eigen::VectorXd free_fluxes_eigen(free_fluxes.length());
    for (int i = 0; i < free_fluxes.length(); ++i) {
        free_fluxes_eigen(i) = free_fluxes[i];
    }

    Matrix all_fluxes_matrix = nullspace * free_fluxes_eigen;


    std::vector<Flux> all_fluxes(reactions.size());
    // fill depended fluxes
    for (int i = 0; i < all_fluxes_matrix.rows(); ++i) {
        all_fluxes[reactions.at(reactions.size() - all_fluxes_matrix.rows() + i).id] = -all_fluxes_matrix(i, 0);
    }

    // fill const fake fluxes
    for (int i = 0; i < reactions.size() - all_fluxes_matrix.rows() - free_fluxes_eigen.size(); ++i) {
        all_fluxes[reactions.at(i).id] = 1;
    }

    // fill free fluxes
    for (int i = 0; i < free_fluxes_eigen.size(); ++i) {
        all_fluxes[reactions.at(reactions.size() - free_fluxes_eigen.size() + i).id] = free_fluxes_eigen[i];
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

Matrix GetRREF(const Matrix &item) {
    Matrix rref_form = item;
    Eigen::FullPivLU<Eigen::MatrixXd> lu_matrix(rref_form);

    for (int pivot = 0; pivot < lu_matrix.rank(); ++pivot) {
        rref_form.row(pivot) /= rref_form(pivot, pivot);
        for (int i = pivot + 1; i < rref_form.rows(); ++i) {
            rref_form.row(i) -= rref_form.row(pivot) * rref_form(i, pivot);
        }
    }

    for (int pivot = lu_matrix.rank() - 1; pivot > 0; --pivot) {
        for (int i = 0; i < pivot; ++i) {
            rref_form.row(i) -= rref_form.row(pivot) * rref_form(i, pivot);
        }
    }

    Matrix A = rref_form.block(0, lu_matrix.rank(), lu_matrix.rank(), rref_form.cols() - lu_matrix.rank());

    return A;
}
