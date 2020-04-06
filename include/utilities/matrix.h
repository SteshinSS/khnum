#pragma once

#include "eigen/Eigen"
#include "eigen/Dense"


namespace khnum {
using Matrix = Eigen::MatrixXd;
using Triplet = Eigen::Triplet<double>;
using SparseMatrix = Eigen::SparseMatrix<double>;
} // namespace khnum