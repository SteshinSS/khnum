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

struct ObjectiveParameters {
    Matrix *nullspace;
    std::vector<EMUNetwork> *networks;
    std::vector<EMUandMID> *input_mids;
    std::vector<EMU> *measured_isotopes;
    std::vector<double> *errors;
    std::vector<EMUandMID> *measurements;
};

std::vector<Flux> EstimateFluxes(const std::vector<EMUNetwork> &networks,
                                           const std::vector<Flux> &initial_fluxes,
                                           const std::vector<FluxVariability> &flux_ranges,
                                           const Matrix &stoichiometry_matrix) {
    // find nullspace
    Eigen::FullPivLU<Matrix> lu_decomposition(stoichiometry_matrix);
    Matrix nullspace = lu_decomposition.kernel();


    // find free fluxes
    const int fluxes_total = stoichiometry_matrix.cols();
    const int nullity = fluxes_total - lu_decomposition.rank();
    alglib::real_1d_array free_fluxes;
    free_fluxes.setlength(nullity);
    for (int i = 0; i < nullity; ++i) {
        free_fluxes[i] = initial_fluxes[fluxes_total - 1 - nullity + i];
    }



    // optimization step
    double epsx = 0.000000001;
    alglib::ae_int_t maxits = 0;
    alglib::minlmstate state;
    alglib::minlmreport report;

    int total_residuals = 0;

}

void CalculateResidual(const alglib::real_1d_array &free_fluxes,
                       alglib::real_1d_array &residuals,
                       void *ptr) {
    ObjectiveParameters *parameters = static_cast<ObjectiveParameters *>(ptr);
    Matrix nullspace = *(parameters->nullspace);

    Eigen::VectorXd free_fluxes_eigen(free_fluxes.length());
    for (int i = 0; i < free_fluxes.length(); ++i) {
        free_fluxes_eigen[i] = free_fluxes[i];
    }

    Matrix all_fluxes = nullspace * free_fluxes_eigen;

    std::vector<Flux> vector_fluxes;
    for (int i = 0; i < all_fluxes.rows(); ++i) {
        vector_fluxes.push_back(all_fluxes(i, 0));
    }



    std::vector<EMUandMID> simulated_mids = CalculateMids(vector_fluxes,
                                                          *(parameters->networks),
                                                          *(parameters->input_mids),
                                                          *(parameters->measured_isotopes));

    int total_residuals = 0;
    for (int isotope = 0; isotope < simulated_mids.size(); ++isotope) {
        for (int mass_shift = 0; mass_shift < simulated_mids[isotope].mid.size(); ++mass_shift) {
            residuals[total_residuals] = simulated_mids[isotope].mid[mass_shift];
            residuals[total_residuals] -= (*(parameters->measurements))[isotope].mid[mass_shift];
            residuals[total_residuals] /= (*(parameters->errors))[total_residuals];
            ++total_residuals;
        }
    }


}


