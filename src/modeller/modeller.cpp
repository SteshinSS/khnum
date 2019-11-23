#include "modeller/modeller.h"

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
#include "modeller/calculate_flux_bounds.h"
#include "modeller/check_model.h"

#include "utilities/debug_utills/debug_prints.h"

namespace khnum {
Modeller::Modeller(ParserResults parser_results) {
    reactions_ = std::move(parser_results.reactions);
    measured_isotopes_ = std::move(parser_results.measured_isotopes);
    measurements_ = std::move(parser_results.measurements);
    excluded_metabolites_ = std::move(parser_results.excluded_metabolites);
    input_substrate_ = std::move(parser_results.input_substrate);
}


void Modeller::CalculateInputSubstrateMids() {
    all_emu_reactions_ = modelling_utills::CreateAllEmuReactions(reactions_, measured_isotopes_);
    input_emu_list_ = modelling_utills::CreateInputEmuList(all_emu_reactions_, input_substrate_);
    input_substrate_mids_ = modelling_utills::CalculateInputMid(input_substrate_, input_emu_list_);
}


void Modeller::CreateEmuNetworks() {
    emu_networks_ = modelling_utills::CreateEmuNetworks(all_emu_reactions_, input_emu_list_, measured_isotopes_);
}


void Modeller::CalculateFluxBounds() {
    modelling_utills::CalculateFluxBounds(reactions_, stoichiometry_matrix_);
    int nullity = nullspace_.cols();

    lower_bounds_.resize(nullity);
    upper_bounds_.resize(nullity);
    for (int i = 0; i < nullity; ++i) {
        lower_bounds_[i] = reactions_[reactions_.size() - nullity + i].computed_lower_bound;
        upper_bounds_[i] = reactions_[reactions_.size() - nullity + i].computed_upper_bound;
    }
}


void Modeller::CreateNullspaceMatrix() {
    std::vector<std::string> full_metabolite_list = modelling_utills::CreateFullMetaboliteList(reactions_);
    std::vector<std::string> included_metabolites = modelling_utills::CreateIncludedMetaboliteList(full_metabolite_list,
                                                                                 excluded_metabolites_);

    stoichiometry_matrix_ = modelling_utills::CreateStoichiometryMatrix(reactions_, included_metabolites);
    nullspace_ = modelling_utills::GetNullspace(stoichiometry_matrix_, reactions_);
    id_to_position_in_depended_fluxes_.resize(reactions_.size());
    const size_t isotopomer_balance_reactions_total = reactions_.size() - nullspace_.rows() - nullspace_.cols();
    for (size_t i = 0; i < reactions_.size(); ++i) {
        if (i < isotopomer_balance_reactions_total) {
            id_to_position_in_depended_fluxes_[reactions_[i].id] = -1;
        } else {
            if (i < isotopomer_balance_reactions_total + nullspace_.rows()) {
                id_to_position_in_depended_fluxes_[reactions_[i].id] = i - isotopomer_balance_reactions_total;
            } else {
                id_to_position_in_depended_fluxes_[reactions_[i].id] = -1;
            }
        }
    }

    for (size_t i = reactions_.size() - nullspace_.cols(); i < reactions_.size(); ++i) {
        free_fluxes_id_.push_back(reactions_[i].id);
    }
}


void Modeller::CalculateMeasurementsCount() {
    measurements_count_ = 0;

    for (const Measurement &measurement : measurements_) {
        measurements_count_ += measurement.mid.size();
    }
}


void Modeller::CheckModelForErrors() {
    modelling_utills::CheckMeasurementsMID(measurements_);
}

Problem Modeller::GetProblem() {
    Problem problem;
    problem.measured_isotopes = measured_isotopes_;
    problem.nullspace = nullspace_;
    problem.measurements = measurements_;
    problem.measurements_count = measurements_count_;
    problem.lower_bounds = lower_bounds_;
    problem.upper_bounds = upper_bounds_;
    problem.reactions_total = reactions_.size();

    GeneratorParameters& simulator_parameters = problem.simulator_parameters_;
    simulator_parameters.networks = emu_networks_;
    simulator_parameters.input_mids = input_substrate_mids_;
    simulator_parameters.measured_isotopes = measured_isotopes_;
    simulator_parameters.nullspace = nullspace_;
    simulator_parameters.free_flux_id_to_nullspace_position = id_to_position_in_depended_fluxes_;
    simulator_parameters.free_fluxes_id = free_fluxes_id_;

    std::vector<ReactionsName> reactions_names;
    for (const Reaction& reaction : reactions_) {
        ReactionsName reaction_name;
        reaction_name.id = reaction.id;
        reaction_name.name = reaction.name;
        reactions_names.push_back(reaction_name);
    }
    problem.reactions = reactions_names;

    return problem;
}
} // namespace khnum