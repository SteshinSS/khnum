#ifndef CFLEX_MODEL_PARSER_H
#define CFLEX_MODEL_PARSER_H

#include "../utilities/reaction_struct.h"

#include <vector>
#include <string>
#include <sstream>
#include <utility>

std::vector<Reaction> ParseReactions(const std::string &model_path);

void FillReaction(Reaction *reaction, std::stringstream &line);

ChemicalEquation ParseChemicalEquation(std::stringstream &line);
Rate ParseRate(const std::string &rate);
ReactionType ParseReactionType(const std::string &type);
std::tuple<Basis, bool> ParseBasis(const std::string &basis);
Deviation ParseDeviation(const std::string &deviation);
Bound ParseLowerBound(const std::string &lower_bound);
Bound ParseUpperBound(const std::string &upper_bound);

ChemicalEquationSide FillEquationSide(const std::string &substrate_equation, const std::string &atom_equation);
ChemicalEquationSide ParseSubstrateEquationSide(const std::string &raw_equation);
void ParseAtomEquationSide(ChemicalEquationSide *equation_side, const std::string &atom_equation);

#endif //CFLEX_MODEL_PARSER_H