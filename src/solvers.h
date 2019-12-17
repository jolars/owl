#pragma once

#include "solvers/solver.h"
#include "solvers/fista.h"
#include "solvers/admm.h"

using namespace arma;

// helper to choose solver
std::unique_ptr<Solver>
setupSolver(const std::string& solver_choice,
            const bool standardize_features,
            const bool is_sparse,
            const bool diagnostics,
            const uword max_passes,
            const uword verbosity,
            const Rcpp::List& args)
{
  using Rcpp::as;
  using namespace arma;

  if (solver_choice == "fista") {

    const double tol_rel_gap          = as<double>(args["tol_rel_gap"]);
    const double tol_infeas           = as<double>(args["tol_infeas"]);

    return std::unique_ptr<FISTA>(new FISTA{standardize_features,
                                            is_sparse,
                                            diagnostics,
                                            max_passes,
                                            tol_rel_gap,
                                            tol_infeas,
                                            verbosity});

  } else {

    const double tol_rel = as<double>(args["tol_rel"]);
    const double tol_abs = as<double>(args["tol_abs"]);
    const double alpha   = as<double>(args["alpha"]);

    return std::unique_ptr<ADMM>(new ADMM{standardize_features,
                                          is_sparse,
                                          diagnostics,
                                          max_passes,
                                          tol_rel,
                                          tol_abs,
                                          alpha,
                                          verbosity});
  }
}
