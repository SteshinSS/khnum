#pragma once

#include "utilities/reaction.h"
#include "utilities/matrix.h"


namespace khnum {
namespace modelling_utills {
Matrix GetNullspace(Matrix matrix, std::vector<Reaction> &reactions);

bool ExchangeRowsToMakePivotNotNull(Matrix &matrix, const int column);

int FindNotNullColumn(const Matrix &matrix, const int currentRow);
} // namespace modelling_utills
} // namespace khnum