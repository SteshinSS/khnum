#pragma once

#include <string>
#include <vector>


using AtomRepresentation = std::string;
using SubstrateCoefficient = double;

struct Substrate {
    std::string name;                 // for example "PYR_EX"
    AtomRepresentation formula;       // "abc"
    SubstrateCoefficient coefficient; // 1.0
};

// need for stl containers
bool operator==(const Substrate &lhs, const Substrate &rhs);


using ChemicalEquationSide = std::vector<Substrate>;

struct ChemicalEquation {
    ChemicalEquationSide left;
    ChemicalEquationSide right;
};


using Rate = double;

/// see FluxPyt User Guide page 3
enum class ReactionType {
    Irreversible,      ///< Irreversible reaction used for both metabolite and isotopomer balance (F)
    Forward,           ///< Forward reaction of reversible one (FR)
    Backward,          ///< Backward reaction of reversible one (R)
    IsotopomerBalance, ///< Reaction used for isotopomer balancing only (S)
    MetaboliteBalance  ///< Irreversible reaction which is not used for only metabolite balance (B)
};

using Basis = double;
using Deviation = double;
using Flux = double;

struct Reaction {
    int id;
    std::string name;
    ChemicalEquation chemical_equation;
    Rate rate;
    ReactionType type;
    Basis basis;
    Deviation deviation;
    bool is_set_free;
    Flux setted_lower_bound;
    Flux setted_upper_bound;
    Flux computed_lower_bound;
    Flux computed_upper_bound;
};