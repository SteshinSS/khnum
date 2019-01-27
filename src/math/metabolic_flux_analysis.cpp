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

};

std::map<std::string, Flux> EstimateFluxes(const std::vector<EMUNetwork> &networks,
                                           const std::map<std::string, Flux> &initial_fluxes,
                                           const std::map<std::string, FluxVariability> &flux_ranges,
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


    // objective function


    // optimization step
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



    std::vector<EMUandMID> simulated_mids = CalculateMids(all_fluxes,
                                                          *(parameters->networks),
                                                          *(parameters->known_mids),
                                                          *(parameters->measured_isotopes));


}


