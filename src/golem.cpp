#include <RcppArmadillo.h>
#include <memory>
#include "solvers.h"
#include "penalties.h"
#include "families.h"

// [[Rcpp::export]]
Rcpp::List
golemDense(arma::mat        X,
           arma::vec        y,
           const Rcpp::List control)
{
  using Rcpp::as;
  using namespace arma;

  // parameter packs for penalty and solver
  auto penalty_args = as<Rcpp::List>(control["penalty_args"]);
  auto solver_args = as<Rcpp::List>(control["solver_args"]);

  auto family_choice = as<std::string>(control["family_choice"]);
  auto standardize = as<std::string>(control["standardize"]);
  auto fit_intercept = as<bool>(control["fit_intercept"]);

  auto family = setupFamily(family_choice);
  auto penalty = setupPenalty(penalty_args);

  // preprocess response using family-specific defaults
  double y_center, y_scale;
  family->preprocessResponse(y, y_center, y_scale, standardize);

  // peprocess featuresa
  rowvec X_center, X_scale;
  preprocessFeatures(X, X_center, X_scale, standardize);

  FISTA solver(solver_args);
  auto result = solver.fit(X, y, family, penalty, fit_intercept);

  auto intercept = as<double>(result["intercept"]);
  auto beta = as<vec>(result["beta"]);
  auto passes = as<uword>(result["passes"]);

  std::tie(intercept, beta) = unstandardize(std::move(intercept),
                                            std::move(beta),
                                            X_center,
                                            X_scale,
                                            y_center,
                                            y_scale,
                                            fit_intercept);

  return Rcpp::List::create(
    Rcpp::Named("intercept")  = intercept,
    Rcpp::Named("beta")       = Rcpp::wrap(beta),
    Rcpp::Named("lambda")     = Rcpp::wrap(penalty->getLambda()),
    Rcpp::Named("sigma")      = penalty->getSigma(),
    Rcpp::Named("passes")     = passes
  );

  return result;
}
