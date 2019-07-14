#pragma once

#include "alglib/ap.h"
#include "eigen/Eigen"

Eigen::VectorXd GetEigenVectorFromAlgLibVector(const alglib::real_1d_array &alglib_vector);
