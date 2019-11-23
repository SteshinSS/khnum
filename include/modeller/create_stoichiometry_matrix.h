#pragma once

#include "utilities/matrix.h"
#include "utilities/reaction.h"


namespace khnum {
namespace modelling_utills {
Matrix CreateStoichiometryMatrix(const std::vector<Reaction> &reactions,
                                 const std::vector<std::string> &metabolite_list);

double GetTotalCoefficient(const ChemicalEquation &chemical_equation, const std::string &metabolite);
} // namespace modelling_utills
} // namespace khnum