#include <RcppArmadillo.h>
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

  auto family_choice = as<std::string>(control["family_choice"]);
  auto penalty_choice = as<std::string>(control["penalty_choice"]);
  auto solver_choice = as<std::string>(control["solver_choice"]);
  auto sigma_type = as<std::string>(control["sigma_type"]);
  auto lambda = as<arma::vec>(control["lambda"]);
  auto lambda_type = as<std::string>(control["lambda_type"]);
  auto sigma = as<double>(control["sigma"]);
  auto standardize = as<std::string>(control["standardize"]);
  auto tol = as<double>(control["tol"]);
  auto fdr = as<double>(control["fdr"]);
  auto max_passes = as<uword>(control["max_passes"]);
  auto fit_intercept = as<bool>(control["fit_intercept"]);

  uword n = X.n_rows;
  uword p = X.n_cols;

  auto family = setupFamily(family_choice, X, y);
  auto penalty = setupPenalty(penalty_choice,
                              lambda,
                              lambda_type,
                              sigma,
                              sigma_type,
                              n,
                              p,
                              fdr);

  // preprocess response using family-specific defaults
  double y_center, y_scale;
  std::tie(y_center, y_scale) = family->preprocessResponse(standardize);

  // peprocess features
  rowvec X_center, X_scale;
  std::tie(X_center, X_scale) = preprocessFeatures(X, standardize);

  FISTA solver(X, y, lambda, family, penalty, fit_intercept, tol, max_passes);

  auto result = solver.fit();

  auto intercept = as<double>(result["intercept"]);
  auto beta = as<vec>(result["beta"]);
  auto passes = as<uword>(result["passes"]);

  std::tie(intercept, beta) = unstandardize(intercept,
                                            beta,
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
