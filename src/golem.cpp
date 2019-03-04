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

  auto n = x.n_rows;
  auto m = y.n_cols;
  auto p = x.n_cols;

  // parameter packs for penalty and solver
  auto penalty_args = as<Rcpp::List>(control["penalty_args"]);
  auto solver_args = as<Rcpp::List>(control["solver_args"]);

  auto family_choice = as<std::string>(control["family_choice"]);
  auto standardize = as<std::string>(control["standardize"]);
  auto fit_intercept = as<bool>(control["fit_intercept"]);
  auto diagnostics = as<bool>(solver_args["diagnostics"]);

  // setup family and response
  auto family = setupFamily(family_choice);
  rowvec y_center(m);
  rowvec y_scale(m);
  family->preprocessResponse(y, y_center, y_scale, standardize);

  // peprocess features
  rowvec x_center(p);
  rowvec x_scale(p);
  preprocessFeatures(x, x_center, x_scale, standardize);

  // setup penalty
  auto penalty = setupPenalty(penalty_args, x, y, y_scale, family);
  auto n_penalties = penalty->pathLength();

  cube betas(p, m, n_penalties);
  cube intercepts(1, m, n_penalties);

  // setup a few return values
  rowvec intercept(m, fill::zeros);
  mat beta(p, m, fill::zeros);

  uvec passes(n_penalties);
  std::vector<std::vector<double>> losses;
  std::vector<std::vector<double>> timings;

  FISTA solver(std::move(intercept),
               std::move(beta),
               solver_args);

  for (uword i = 0; i < n_penalties; ++i) {
    Results res = solver.fit(x, y, family, penalty, fit_intercept);

    betas.slice(i) = res.beta;
    intercepts.slice(i) = res.intercept;
    passes(i) = res.passes;

    if (diagnostics) {
      losses.push_back(res.loss);
      timings.push_back(res.time);
    }
  }

  unstandardize(intercepts,
                betas,
                x_center,
                x_scale,
                y_center,
                y_scale,
                fit_intercept);

  return Rcpp::List::create(
    Named("intercept")   = wrap(intercepts),
    Named("beta")        = wrap(betas),
    Named("penalty")     = wrap(penalty->getParams(y_scale)),
    Named("passes")      = passes,
    Named("loss")        = wrap(losses),
    Named("time")        = wrap(timings)
  );
}
