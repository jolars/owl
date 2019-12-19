#pragma once

#include "solvers/solver.h"
#include "solvers/fista.h"

using namespace arma;
using namespace Rcpp;

// helper to choose solver
std::unique_ptr<Solver> setupSolver(const std::string& solver_choice,
                                    const bool standardize_features,
                                    const bool is_sparse,
                                    const bool diagnostics,
                                    const uword max_passes,
                                    const double tol_rel_gap,
                                    const double tol_infeas,
                                    const uword verbosity)
{
  return std::unique_ptr<FISTA>(new FISTA{standardize_features,
                                          is_sparse,
                                          diagnostics,
                                          max_passes,
                                          tol_rel_gap,
                                          tol_infeas,
                                          verbosity});
}
