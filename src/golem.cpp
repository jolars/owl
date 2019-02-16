#include <RcppArmadillo.h>
#include "solvers.h"
#include "penalties.h"
#include "families.h"

// [[Rcpp::export]]
Rcpp::List
denseGolem(const arma::mat& X,
           const arma::vec& y,
           const Rcpp::List control)
{
  using Rcpp::as;

  auto family_choice = as<std::string>(control["family_choice"]);
  auto penalty_choice = as<std::string>(control["penalty_choice"]);
  auto solver_choice = as<std::string>(control["solver_choice"]);
  auto lambda = as<arma::vec>(control["lambda"]);

  auto family = setupFamily(family_choice, X, y);
  auto penalty = setupPenalty(penalty_choice, lambda);

  FISTA solver(X, y, lambda, family, penalty);

  auto result = solver.fit();

  return result;
}
