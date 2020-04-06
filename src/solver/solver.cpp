#include "solver/solver.h"

#include <random>
#include <chrono>
#include <ctime>
#include "alglib/optimization.h"

#include "simulator/simulator.h"
#include "utilities/debug_utills/debug_prints.h"
#include "utilities/get_eigen_vec_from_alglib_vec.h"
#include "simulator/generator.h"



namespace khnum {


Solver::Solver(const Problem &problem, const SimulatorGenerator &generator) {
    new_simulator_.emplace(generator.Generate());
    reactions_num_ = problem.reactions_total;
    nullspace_ = problem.nullspace;
    measured_mids_ = problem.measurements;
    measurements_count_ = problem.measurements_count;
    use_analytic_gradient_ = problem.use_analytic_jacobian;
    reactions_ = problem.reactions;

    nullity_ = nullspace_.cols();
    free_fluxes_.setlength(nullity_);

    iteration_ = 0;
    iteration_total_ = 10;

    lower_bounds_.setlength(nullity_);
    upper_bounds_.setlength(nullity_);
    for (int i = 0; i < nullity_; ++i) {
        lower_bounds_[i] = problem.lower_bounds[i];
        upper_bounds_[i] = problem.upper_bounds[i];
    }

}


std::vector<alglib::real_1d_array> Solver::GetResult() {
    return all_solutions_;
}


void Solver::Solve() {
    std::random_device randomizer;
    std::mt19937 random_source(randomizer());

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    for (iteration_ = 0; iteration_ < iteration_total_; ++iteration_) {
        GenerateInitialPoints(random_source);
        if (iteration_ == 0) {
            SetOptimizationParameters();
        } else {
            alglib::minlmrestartfrom(state_, free_fluxes_);
        }

        alglib::real_1d_array new_solution = RunOptimization();

        all_solutions_.emplace_back(new_solution);
    }
    end = std::chrono::system_clock::now();
    double elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>
        (end-start).count();

    elapsed_milliseconds /= 1000;

    std::cout << "Average time: " << static_cast<double>(elapsed_milliseconds) / iteration_total_
       << " seconds per iteration" << std::endl;
}


void Solver::GenerateInitialPoints(std::mt19937 &random_source) {
    std::uniform_real_distribution<> get_random_point(0.0, 1.0);

    for (int i = 0; i < nullity_; ++i) {
        free_fluxes_[i] = lower_bounds_[i] + get_random_point(random_source) * (upper_bounds_[i] - lower_bounds_[i]);
    }

    PrintStartMessage();
}


void Solver::SetOptimizationParameters() {
    alglib::ae_int_t maxits = 500;
    const double epsx = 0.000001;

    if (use_analytic_gradient_) {
        alglib::minlmcreatevj(nullity_, measurements_count_, free_fluxes_, state_);
        alglib::minlmoptguardgradient(state_, 0.000);
    } else {
        alglib::minlmcreatev(nullity_, measurements_count_, free_fluxes_, 0.001, state_);
    }

    alglib::minlmsetacctype(state_, 1);
    alglib::minlmsetcond(state_, epsx, maxits);
    alglib::minlmsetbc(state_, lower_bounds_, upper_bounds_);

    SetConstraints();
}


// Set constraints that Vdep = -nullspace * Vfree > 0
// This is the same as nullspace * VFree < 0
void Solver::SetConstraints() {
    alglib::real_2d_array constraint;
    constraint.setlength(nullspace_.rows(), nullspace_.cols() + 1);
    for (int row = 0; row < nullspace_.rows(); ++row) {
        for (int col = 0; col < nullspace_.cols(); ++col) {
            constraint(row, col) = nullspace_(row, col);
        }
        constraint(row, nullspace_.cols()) = 0;
    }

    alglib::integer_1d_array types;
    types.setlength(nullspace_.rows());
    for (int i = 0; i < nullspace_.rows(); ++i) {
        types[i] = -1;
    }

    alglib::minlmsetlc(state_, constraint, types);
}


void Solver::PrintStartMessage() {
    //std::cout << "Start " << iteration_ << " iteration from: " << std::endl;
    /*
    for (int i = 0; i < nullity_; ++i) {
        std::cout << reactions_[reactions_num_ - nullity_ + i].id + 1 << ": " <<  free_fluxes_[i] << std::endl;
    }
    std::cout << std::endl; */
}


alglib::real_1d_array Solver::RunOptimization() {
    if (use_analytic_gradient_) {
        alglib::minlmoptimize(state_, AlglibCallback, JacobianCallback, nullptr, this, alglib::xdefault);

        alglib::optguardreport ogrep;
        alglib::minlmoptguardresults(state_, ogrep);
        if (ogrep.badgradsuspected) {
            std::cout << "Bad jacobian. Should be: " << std::endl;
            for (int i = 0; i < ogrep.badgradnum.rows(); ++i) {
                for (int j = 0; j < ogrep.badgradnum.cols(); ++j) {
                    std::cout << ogrep.badgradnum(i, j) << " ";
                }
                std::cout << std::endl;
            }

            std::cout << std::endl << "But actually: " << std::endl;
            for (int i = 0; i < ogrep.badgraduser.rows(); ++i) {
                for (int j = 0; j < ogrep.badgraduser.cols(); ++j) {
                    std::cout << ogrep.badgraduser(i, j) << " ";
                }
                std::cout << std::endl;
            }

            std::cout << "Fluxes: " << std::endl;
            for (int i = 0; i < ogrep.badgradxbase.length(); ++i) {
                std::cout << ogrep.badgradxbase(i) << "  ";
            }
            std::cout << std::endl;
        }
    } else {
        alglib::minlmoptimize(state_, AlglibCallback, nullptr, this, alglib::xdefault);
    }

    alglib::real_1d_array final_free_fluxes;
    alglib::minlmresults(state_, final_free_fluxes, report_);

    PrintFinalMessage(final_free_fluxes);

    return final_free_fluxes;
}

void JacobianCallback(const alglib::real_1d_array &free_fluxes,
                      alglib::real_1d_array &fi,
                      alglib::real_2d_array &jac, void *ptr) {
    Solver* solver = static_cast<Solver*>(ptr);
    solver->in_jacobian = true;
    AlglibCallback(free_fluxes, fi, ptr);
    solver->in_jacobian = false;
    solver->FillJacobian(jac);
}


void AlglibCallback(const alglib::real_1d_array &free_fluxes,
                    alglib::real_1d_array &residuals, void *ptr) {
    Solver* solver = static_cast<Solver*>(ptr);
    solver->CalculateResidual(free_fluxes, residuals);
}



void Solver::FillJacobian(alglib::real_2d_array &jac) {
    int total_residuals = 0;
    for (size_t isotope = 0; isotope < measured_mids_.size(); ++isotope) {
        for (size_t mass_shift = 0; mass_shift < measured_mids_[isotope].errors.size(); ++mass_shift) {
            for (size_t flux = 0; flux < diff_results_.size(); ++flux) {
                jac(total_residuals, flux) = diff_results_[flux][isotope].mid[mass_shift] / measured_mids_[isotope].errors[mass_shift];
            }
            ++total_residuals;
        }
    }
}


void Solver::CalculateResidual(const alglib::real_1d_array &free_fluxes,
                               alglib::real_1d_array &residuals) {
    std::vector<Flux> calculated_fluxes = CalculateAllFluxesFromFree(free_fluxes);
    SimulatorResult result = new_simulator_->CalculateMids(calculated_fluxes, in_jacobian);

    diff_results_ = result.diff_results;

    Fillf0Array(residuals, result.simulated_mids);
}


std::vector<Flux> Solver::CalculateAllFluxesFromFree(const alglib::real_1d_array &free_fluxes_alglib) {
    Eigen::VectorXd free_fluxes = GetEigenVectorFromAlgLibVector(free_fluxes_alglib);
    Matrix depended_fluxes_matrix = -nullspace_ * free_fluxes;
    std::vector<Flux> all_fluxes(reactions_num_, -1);
    // non metabolite balance reactions
    const int depended_reactions_total = depended_fluxes_matrix.rows();
    const int
        isotopomer_balance_reactions_total = reactions_num_ - depended_reactions_total - free_fluxes_alglib.length();

    for (int i = 0; i < isotopomer_balance_reactions_total; ++i) {
        all_fluxes[reactions_.at(i).id] = 1;
    }

    for (int i = 0; i < depended_reactions_total; ++i) {
        all_fluxes[reactions_.at(i + isotopomer_balance_reactions_total).id] = depended_fluxes_matrix(i, 0);
    }

    for (int i = 0; i < free_fluxes_alglib.length(); ++i) {
        all_fluxes[reactions_.at(reactions_num_ - free_fluxes_alglib.length() + i).id] = free_fluxes[i];
    }

    return all_fluxes;
}


void Solver::Fillf0Array(alglib::real_1d_array &residuals, const std::vector<EmuAndMid> &simulated_mids) {
    int total_residuals = 0;
    for (size_t isotope = 0; isotope < simulated_mids.size(); ++isotope) {
        for (size_t mass_shift = 0; mass_shift < simulated_mids[isotope].mid.size(); ++mass_shift) {
            residuals[total_residuals] = simulated_mids[isotope].mid[mass_shift];
            residuals[total_residuals] -= (measured_mids_[isotope].mid[mass_shift]);
            residuals[total_residuals] /= measured_mids_[isotope].errors[mass_shift];
            ++total_residuals;
        }
    }
}


double Solver::GetSSR(const alglib::real_1d_array &residuals) {
    double answer = 0.0;
    for (int measurement = 0; measurement < measurements_count_; ++measurement) {
        answer += residuals(measurement) * residuals(measurement);
    }
    return answer;
}


void Solver::PrintFinalMessage(const alglib::real_1d_array &free_fluxes) {
    std::vector<Flux> final_all_fluxes = CalculateAllFluxesFromFree(free_fluxes);
    SimulatorResult result = new_simulator_->CalculateMids(final_all_fluxes, in_jacobian);
    alglib::real_1d_array residuals;
    residuals.setlength(measurements_count_);
    Fillf0Array(residuals, result.simulated_mids);

    double ssr = GetSSR(residuals);
/*
    std::cout << "Finish at: " << std::endl;

    for (int i = 0; i < free_fluxes.length(); ++i) {
        std::cout << reactions_[reactions_num_ - free_fluxes.length() + i].id + 1 <<
                  " = " << free_fluxes[i] << std::endl;
    } */
    std::cout << " Finish with SSR: " << ssr << " in " << report_.iterationscount << " steps." << std::endl;
}
} // namespace khnum
