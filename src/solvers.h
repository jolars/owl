#pragma once

#include "solvers/solver.h"
#include "solvers/fista.h"

// helper to choose solver
std::unique_ptr<Solver>
setupSolver(const std::string& solver_choice,
            const bool standardize_features,
            const bool is_sparse,
            const Rcpp::List& args)
{
  return std::unique_ptr<FISTA>(new FISTA{standardize_features,
                                          is_sparse,
                                          args});
}
