#include "modeller/create_nullspace.h"

#include "utilities/reaction.h"
#include "utilities/matrix.h"

#include "iostream"


namespace khnum {
namespace modelling_utills {
const double epsilon = 1.e-12;
// Transform stoichiometry matrix to form
/*
 *
 * |1 0 0 0 0 a1 a2 |
 * |0 1 0 0 0 a3 a4 |
 * |0 0 1 0 0 a5 a6 | = | I  A |
 * |0 0 0 1 0 a7 a8 |
 * |0 0 0 0 1 a9 a10|
 */
// So Vdep = -A * Vfree
// And return the A matrix

Matrix GetNullspace(Matrix matrix, std::vector<Reaction> &reactions) {
    // As we work with the diagonal elements, pivot.column == pivot.row
    // So below I use it as synonymous


    const int metabolite_balance_reactions_total = reactions.size() - matrix.cols();

    for (int column = 0; column < matrix.rows(); ++column) {
        if (abs(matrix(column, column)) < epsilon) {
            // Try to exchange with row below
            bool is_found_not_null_pivot = ExchangeRowsToMakePivotNotNull(matrix, column);

            if (!is_found_not_null_pivot) {
                // We need to exchange columns, so as the fluxes order
                const int column_to_swap = FindNotNullColumn(matrix, column);
                if (column_to_swap == -1) {
                    throw std::runtime_error(
                        "Can't transform stoichiometry matrix at column number " + std::to_string(column));
                }
                matrix.col(column).swap(matrix.col(column_to_swap));
                std::swap(reactions[metabolite_balance_reactions_total + column],
                          reactions[metabolite_balance_reactions_total + column_to_swap]);
                std::cout << "Reaction num " << metabolite_balance_reactions_total + column << " and num "
                          << metabolite_balance_reactions_total + column_to_swap << " has swapped" << std::endl;
                if (matrix(column, column) == 0) {
                    ExchangeRowsToMakePivotNotNull(matrix, column);
                }
            }
        }

        // Now it's guaranteed pivot is not null
        matrix.row(column) /= matrix(column, column);
        for (int row = column + 1; row < matrix.rows(); ++row) {
            matrix.row(row) -= matrix.row(column) * matrix(row, column);
        }

    }

    // Now we have upper-diagonal matrix with units at diagonal
    for (int column = matrix.rows() - 1; column >= 0; --column) {
        for (int row = column - 1; row >= 0; --row) {
            matrix.row(row) -= matrix.row(column) * matrix(row, column);
        }
    }

    Matrix result = matrix.block(0, matrix.rows(), matrix.rows(), matrix.cols() - matrix.rows());
    return result;
}


bool ExchangeRowsToMakePivotNotNull(Matrix &matrix, const int column) {
    for (int row = column + 1; row < matrix.rows(); ++row) {
        if (abs(matrix(row, column)) > epsilon) {
            matrix.row(column).swap(matrix.row(row));
            return true;
        }
    }
    return false;
}


int FindNotNullColumn(const Matrix &matrix, const int currentRow) {
    for (int column = matrix.rows() - 1; column < matrix.cols(); ++column) {
        for (int row = currentRow; row < matrix.rows(); ++row) {
            if (abs(matrix(row, column)) > epsilon) {
                return column;
            }
        }
    }
    return -1;
}
} // namespace modelling_utills
} // namespace khnum