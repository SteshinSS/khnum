#pragma once

#include <string>
#include <vector>


namespace khnum {
using Fraction = double;

struct Mixture {
    std::vector<Fraction> fractions;
    double ratio;
};

struct InputSubstrate {
    std::string name;
    std::vector<Mixture> mixtures;
};
} // namespace khnum