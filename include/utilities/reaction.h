#pragma once

#include <string>
#include <optional>
#include <vector>


namespace khnum {
using AtomRepresentation = std::string;
using SubstrateCoefficient = double;

struct Substrate {
    int size; // 3
    int id; // 0
    std::string name;                 // for example, "PYR_EX"
    AtomRepresentation formula;       // "abc"
    SubstrateCoefficient substrate_coefficient_ = 0.0; // 1.0
};

// need for stl containers
bool operator==(const Substrate &lhs, const Substrate &rhs);

using ChemicalEquationSide = std::vector<Substrate>;

struct AtomTransition {
  int substrate_pos;
  int product_pos;
  int substrate_atom;
  int product_atom;
};

struct ChemicalEquation {
    ChemicalEquationSide left;
    ChemicalEquationSide right;
    std::vector<AtomTransition> atom_transitions;
};




/// see OpenFlux paper
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

    ReactionType type;
    Basis basis;
    Deviation deviation;
    bool is_set_free;
    std::optional<Flux> setted_lower_bound;
    std::optional<Flux> setted_upper_bound;
    Flux computed_lower_bound;
    Flux computed_upper_bound;
};
} // namespace khnum