#ifndef CFLEX_INPUT_SUBSTRATE_H
#define CFLEX_INPUT_SUBSTRATE_H

#include <string>
#include <vector>

using Fraction = double;

struct Mixture {
  std::vector<Fraction> fractions;
  double ratio;
};

struct InputSubstrate {
  std::string name;
  std::vector<Mixture> mixtures;
};

#endif //CFLEX_INPUT_SUBSTRATE_H
