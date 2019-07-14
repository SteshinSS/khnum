#pragma once

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

