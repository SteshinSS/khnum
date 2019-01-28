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


std::map<std::string, Flux> EstimateFluxes(ObjectiveParameters *parameters,
                                           const std::map<std::string, Flux> &initial_fluxes,
                                           const std::map<std::string, FluxVariability> &flux_ranges,
                                           const Matrix &stoichiometry_matrix,
                                           const std::vector<Reaction> &reactions) {
    // find nullspace
    Eigen::FullPivLU<Matrix> lu_decomposition(stoichiometry_matrix);
    Matrix nullspace = lu_decomposition.kernel();
    parameters->nullspace = &nullspace;

    // find free fluxes
    const int fluxes_total = stoichiometry_matrix.cols();
    const int nullity = fluxes_total - lu_decomposition.rank();
    alglib::real_1d_array free_fluxes;
    alglib::real_1d_array lower_bounds;
    alglib::real_1d_array upper_bounds;

    free_fluxes.setlength(nullity);
    lower_bounds.setlength(nullity);
    upper_bounds.setlength(nullity);

    const int reaction_total = reactions.size();
    for (int i = 0; i < nullity; ++i) {
        free_fluxes[i] = initial_fluxes.at(reactions[reaction_total - 1 - nullity + i].name);
        lower_bounds[i] = (flux_ranges.at(reactions[reaction_total  - 1 - nullity + i].name)).lower_bound;
        upper_bounds[i] = (flux_ranges.at(reactions[reaction_total - 1 - nullity + i].name)).upper_bound;
    }



    // optimization step
    double epsx = 0.000000001;
    alglib::ae_int_t maxits = 0;
    alglib::minlmstate state;
    alglib::minlmreport report;

    int measurements_count = 0;
    for (const Measurement &measurement : *(parameters->measurements)) {
        measurements_count += measurement.mid.size();
    }

    alglib::minlmcreatev(measurements_count, free_fluxes, 0.00001, state);
    alglib::minlmsetcond(state, epsx, maxits);
    alglib::minlmsetbc(state, lower_bounds, upper_bounds);

    alglib::minlmoptimize(state, CalculateResidual, NULL, parameters);
    alglib::minlmresults(state, free_fluxes, report);

    Eigen::VectorXd free_fluxes_eigen(free_fluxes.length());
    for (int i = 0; i < free_fluxes.length(); ++i) {
        free_fluxes_eigen[i] = free_fluxes[i];
    }

    Matrix all_fluxes = nullspace * free_fluxes_eigen;

    std::map<std::string, Flux> real_fluxes;
    for (int i = 0; i < all_fluxes.rows(); ++i) {
        real_fluxes[(parameters->reactions)->at(i).name] = all_fluxes(i, 0);
    }

    for (int i = all_fluxes.rows(); i < reaction_total; ++i) {
        real_fluxes[(parameters->reactions)->at(i).name] = free_fluxes_eigen(i - all_fluxes.rows());
    }

    return real_fluxes;
}

void CalculateResidual(const alglib::real_1d_array &free_fluxes,
                       alglib::real_1d_array &residuals,
                       void *ptr) {
    ObjectiveParameters &parameters = *(static_cast<ObjectiveParameters *>(ptr));
    Matrix nullspace = *(parameters.nullspace);

    Eigen::VectorXd free_fluxes_eigen(free_fluxes.length());
    for (int i = 0; i < free_fluxes.length(); ++i) {
        free_fluxes_eigen[i] = free_fluxes[i];
    }

    Matrix all_fluxes = nullspace * free_fluxes_eigen;

    std::map<std::string, Flux> real_fluxes;
    for (int i = 0; i < all_fluxes.rows(); ++i) {
        real_fluxes[(parameters.reactions)->at(i).name] = all_fluxes(i, 0);
    }
    for (int i = all_fluxes.rows(); i < parameters.reactions->size(); ++i) {
        real_fluxes[(parameters.reactions)->at(i).name] = free_fluxes_eigen(i - all_fluxes.rows());
    }


    std::vector<EMUandMID> simulated_mids = CalculateMids(real_fluxes,
                                                          *(parameters.networks),
                                                          *(parameters.input_mids),
                                                          *(parameters.measured_isotopes));

    int total_residuals = 0;
    for (int isotope = 0; isotope < simulated_mids.size(); ++isotope) {
        for (int mass_shift = 0; mass_shift < simulated_mids[isotope].mid.size(); ++mass_shift) {
            residuals[total_residuals] = simulated_mids[isotope].mid[mass_shift];
            residuals[total_residuals] -= (*(parameters.measurements))[isotope].mid[mass_shift];
            residuals[total_residuals] /= (*(parameters.measurements))[isotope].errors[mass_shift];
            ++total_residuals;
        }
    }

    return;
}


