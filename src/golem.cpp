#include <RcppArmadillo.h>
#include <memory>
#include "solvers.h"
#include "penalties.h"
#include "families.h"

// [[Rcpp::export]]
Rcpp::List
golemDense(arma::mat x,
           arma::mat y,
           const Rcpp::List control)
{
  using namespace arma;

  using Rcpp::as;
  using Rcpp::Named;
  using Rcpp::wrap;

  // auto n = x.n_rows;
  auto m = y.n_cols;
  auto p = x.n_cols;

  // parameter packs for penalty and solver
  auto penalty_args = as<Rcpp::S4>(control["penalty"]);
  auto solver_args = as<Rcpp::S4>(control["solver"]);
  auto lipschitz_constant = as<double>(control["lipschitz_constant"]);

  auto family_args = as<Rcpp::S4>(control["family"]);
  auto fit_intercept = as<bool>(control["fit_intercept"]);
  bool diagnostics = solver_args.slot("diagnostics");

  // // setup family and response
  auto family = setupFamily(family_args.slot("name"));
  auto penalty = setupPenalty(penalty_args);
  auto n_penalties = penalty->pathLength();

  cube betas(p, m, n_penalties);
  cube intercepts(1, m, n_penalties);

  // setup a few return values
  rowvec intercept(m, fill::zeros);
  mat beta(p, m, fill::zeros);

  uvec passes(n_penalties);
  std::vector<std::vector<double>> primals;
  std::vector<std::vector<double>> duals;
  std::vector<std::vector<double>> timings;

  FISTA solver(std::move(intercept),
               std::move(beta),
               lipschitz_constant,
               solver_args);

  for (uword i = 0; i < n_penalties; ++i) {
    Results res = solver.fit(x, y, family, penalty, fit_intercept);

    betas.slice(i) = res.beta;
    intercepts.slice(i) = res.intercept;
    passes(i) = res.passes;
    penalty->step(i + 1); // selects next lambda for LASSO for instance

    if (diagnostics) {
      primals.push_back(res.primals);
      duals.push_back(res.duals);
      timings.push_back(res.time);
    }
  }

  return Rcpp::List::create(
    Named("intercept")   = wrap(intercepts),
    Named("beta")        = wrap(betas),
    Named("passes")      = passes,
    Named("primals")     = wrap(primals),
    Named("duals")       = wrap(duals),
    Named("time")        = wrap(timings)
  );
}
