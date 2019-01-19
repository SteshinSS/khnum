#ifndef CFLEX_REACTION_STRUCT_H_
#define CFLEX_REACTION_STRUCT_H_

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
using Bound = double;



struct Reaction {
  std::string name;
  ChemicalEquation chemical_equation;
  Rate rate;
  ReactionType type;
  Basis basis;
  Deviation deviation;
  Bound lower_bound;
  Bound upper_bound;
  bool is_prime_basis;
  bool is_basis_x;
};

// probably better to move into a separate header
using Flux = double;

#endif