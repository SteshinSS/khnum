#include "modeller.h"

#include <math.h>
#include <iostream>

#include "utilities/matrix.h"
#include "modeller/create_emu_networks.h"
#include "modeller/calculate_input_mid.h"
#include "modeller/create_emu_reactions.h"
#include "modeller/create_emu_list.h"
#include "modeller/create_metabolite_list.h"
#include "modeller/create_stoichiometry_matrix.h"
#include "modeller/create_nullspace.h"

Modeller::Modeller(ParserResults parser_results) {
    reactions_ = std::move(parser_results.reactions);
    measured_isotopes_ = std::move(parser_results.measured_isotopes);
    measurements_ = std::move(parser_results.measurements);
    excluded_metabolites_ = std::move(parser_results.excluded_metabolites);
    input_substrate_ = std::move(parser_results.input_substrate);
}

void Modeller::CalculateInputSubstrateMids() {
    all_emu_reactions_ = CreateAllEMUReactions(reactions_, measured_isotopes_);
    input_emu_list_ = CreateInputEMUList(all_emu_reactions_, input_substrate_);
    input_substrate_mids_ = CalculateInputMid(input_substrate_, input_emu_list_);
}

void Modeller::CreateEmuNetworks() {
    emu_networks_ = CreateEMUNetworks(all_emu_reactions_, input_emu_list_, measured_isotopes_);
}

void Modeller::CalculateFluxBounds() {
    for (Reaction &reaction : reactions_) {
        if (std::isnan(reaction.basis)) {
            reaction.computed_upper_bound = reaction.setted_upper_bound;
            reaction.computed_lower_bound = reaction.setted_lower_bound;
        } else {
            if (std::isnan(reaction.deviation)) {
                reaction.computed_upper_bound = reaction.basis;
                reaction.computed_lower_bound = reaction.basis;
            } else {
                reaction.computed_upper_bound = reaction.basis + reaction.deviation;
                reaction.computed_lower_bound = reaction.basis - reaction.deviation;
            }
        }
    }
}

void Modeller::CreateNullspaceMatrix() {
    std::vector<std::string> full_metabolite_list = CreateFullMetaboliteList(reactions_);
    std::vector<std::string> included_metabolites = CreateIncludedMetaboliteList(full_metabolite_list,
                                                                                 excluded_metabolites_);

    Matrix stoichiometry_matrix = CreateStoichiometryMatrix(reactions_, included_metabolites);
    std::cout << stoichiometry_matrix << std::endl << std::endl;

    nullspace_ = GetNullspace(stoichiometry_matrix, reactions_);

    std::cout << nullspace_ << std::endl << std::endl;
}


void Modeller::CalculateMeasurementsCount() {
    measurements_count_ = 0;

    for (const Measurement &measurement : measurements_) {
        measurements_count_ += measurement.mid.size();
    }
}


Problem Modeller::GetProblem() {
    Problem problem;
    problem.reactions = reactions_;
    problem.measured_isotopes = measured_isotopes_;
    problem.nullspace = nullspace_;
    problem.networks = emu_networks_;
    problem.input_mids = input_substrate_mids_;
    problem.measurements = measurements_;
    problem.measurements_count = measurements_count_;

    return problem;
}


