#include "modeller/create_nullspace.h"

#include "utilities/reaction.h"
#include "utilities/matrix.h"

#include "iostream"

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

Matrix GetNullspace(Matrix matrix, std::vector<Reaction>& reactions) {
    // As we work with the diagonal elements, pivot.column == pivot.row
    // So below I use it as synonymous
    const int metabolite_balance_reactions_total = reactions.size() - matrix.cols();


    for (int column = 0; column < matrix.rows(); ++column) {
        if (matrix(column, column) == 0) {
            // Try to exchange with row below
            bool isFoundNotNullPivot = ExchangeRowsToMakePivotNotNull(matrix, column);

            if (!isFoundNotNullPivot) {
                // We need to exchange columns, so as the fluxes order
                const int columnToSwap = FindNotNullColumn(matrix, column);
                if (columnToSwap == -1) {
                    throw std::runtime_error("Can't transform stoichiometry matrix at column number " + std::to_string(column));
                }
                matrix.col(column).swap(matrix.col(columnToSwap));
                std::swap(reactions[metabolite_balance_reactions_total + column], reactions[metabolite_balance_reactions_total + columnToSwap]);
                std::cout << "Reaction num " << metabolite_balance_reactions_total + column << " and num " << metabolite_balance_reactions_total + columnToSwap << " has swapped" << std::endl;
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

bool ExchangeRowsToMakePivotNotNull(Matrix& matrix, const int column) {
    for (int row = column + 1; row < matrix.rows(); ++row) {
        if (matrix(row, column) != 0) {
            matrix.row(column).swap(matrix.row(row));
            return true;
        }
    }
    return false;
}

int FindNotNullColumn(const Matrix& matrix, const int currentRow) {
    for (int column = matrix.rows() - 1; column < matrix.cols(); ++column) {
        for (int row = currentRow; row < matrix.rows(); ++row) {
            if (matrix(row, column) != 0) {
                return column;
            }
        }
    }
    return -1;
}