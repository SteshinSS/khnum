#include "utilities/get_eigen_vec_from_alglib_vec.h"

#include "alglib/ap.h"
#include "eigen/Eigen"


namespace khnum {
Eigen::VectorXd GetEigenVectorFromAlgLibVector(const alglib::real_1d_array &alglib_vector) {
    Eigen::VectorXd eigen_vector(alglib_vector.length());
    for (int i = 0; i < alglib_vector.length(); ++i) {
        eigen_vector(i) = alglib_vector[i];
    }

    return eigen_vector;
}
} //namespace khnum