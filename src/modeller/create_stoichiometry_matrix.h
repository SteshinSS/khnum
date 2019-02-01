#ifndef CFLEX_CREATE_STOICHIOMETRY_MATRIX_H
#define CFLEX_CREATE_STOICHIOMETRY_MATRIX_H

#include "math_utilites.h"
#include "reaction_struct.h"

Matrix CreateStoichiometryMatrix(const std::vector<Reaction> &reactions,
                                 const std::vector<std::string> &metabolite_list);

double GetTotalCoefficient(const ChemicalEquation &chemical_equation, const std::string &metabolite);
#endif //CFLEX_CREATE_STOICHIOMETRY_MATRIX_H
