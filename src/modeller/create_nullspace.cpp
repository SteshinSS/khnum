#include "create_nullspace.h"

#include "../math/math_utilites.h"

Matrix GetNullspace(const Matrix &matrix) {
    Matrix nullspace = matrix.fullPivLu().kernel();
    nullspace = GetRREF(nullspace);
    return nullspace;
}

Matrix GetRREF(const Matrix &item) {
    Matrix result = item;
    for (int pivot = item.cols() - 1; pivot > 0; --pivot) {
        result.col(pivot) /= result(item.rows() - item.cols() + pivot, pivot);
        for (int i = 0; i < pivot; ++i) {
            result.col(i) -= result.col(pivot) * result(item.rows() - item.cols() + pivot, i);
        }
    }
    result.col(0) /= result(item.rows() - item.cols(), 0);
    for (int pivot = 0; pivot < item.cols(); ++pivot) {
        for (int i = pivot + 1; i < item.cols(); ++i) {
            result.col(i) -= result.col(pivot) * result(item.rows() - item.cols() + pivot, i);
        }
    }

    return result;
}
