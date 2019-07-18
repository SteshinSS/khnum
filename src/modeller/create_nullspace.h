#pragma once

#include "utilities/reaction.h"
#include "utilities/matrix.h"

Matrix GetNullspace(Matrix matrix, std::vector<Reaction>& reactions);

bool ExchangeRowsToMakePivotNotNull(Matrix& matrix, const int column);

int FindNotNullColumn(const Matrix& matrix, const int currentRow);
