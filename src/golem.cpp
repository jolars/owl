#include <RcppArmadillo.h>
#include <memory>
#include "solvers.h"
#include "penalties.h"
#include "families.h"

template <typename T>
Rcpp::List
golemCpp(const T& x,
         const arma::mat& y,
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
  auto penalty_args = as<Rcpp::List>(control["penalty"]);
  auto solver_args = as<Rcpp::List>(control["solver"]);
  auto groups = as<Rcpp::List>(control["groups"]);
  auto lipschitz_constant = as<double>(control["lipschitz_constant"]);

  auto family_args = as<Rcpp::List>(control["family"]);
  auto fit_intercept = as<bool>(control["fit_intercept"]);
  auto diagnostics = as<bool>(solver_args["diagnostics"]);
  auto standardize_features = as<bool>(control["standardize_features"]);
  auto is_sparse = as<bool>(control["is_sparse"]);
  auto n_penalties = as<uword>(control["n_penalties"]);

  // get scaled vector of feature matrix centers for use in sparse fitting
  vec x_scaled_center = as<vec>(control["x_scaled_center"]);

  // setup family and response
  auto family = setupFamily(as<std::string>(family_args["name"]),
                            fit_intercept,
                            standardize_features);

  auto penalty = setupPenalty(penalty_args, groups);

  cube betas(p, m, n_penalties);
  cube intercepts(1, m, n_penalties);

  // initialize estimates
  auto intercept_init = as<rowvec>(control["intercept_init"]);
  auto beta_init      = as<mat>(control["beta_init"]);

  uvec passes(n_penalties);
  std::vector<std::vector<double>> primals;
  std::vector<std::vector<double>> duals;
  std::vector<std::vector<double>> timings;
  std::vector<std::vector<double>> infeasibilities;

  FISTA solver(intercept_init,
               beta_init,
               lipschitz_constant,
               standardize_features,
               x_scaled_center,
               is_sparse,
               solver_args);

  for (uword path_iter = 0; path_iter < n_penalties; ++path_iter) {
    Results res = solver.fit(x, y, family, penalty, fit_intercept, path_iter);

    betas.slice(path_iter) = res.beta;
    intercepts.slice(path_iter) = res.intercept;
    passes(path_iter) = res.passes;

    if (diagnostics) {
      primals.push_back(res.primals);
      duals.push_back(res.duals);
      infeasibilities.push_back(res.infeasibilities);
      timings.push_back(res.time);
    }
  }

  return Rcpp::List::create(
    Named("intercept")       = wrap(intercepts),
    Named("beta")            = wrap(betas),
    Named("passes")          = passes,
    Named("primals")         = wrap(primals),
    Named("duals")           = wrap(duals),
    Named("infeasibilities") = wrap(infeasibilities),
    Named("time")            = wrap(timings)
  );
}


// [[Rcpp::export]]
Rcpp::List
golemSparse(const arma::sp_mat& x,
            const arma::mat& y,
            const Rcpp::List control)
{
  return golemCpp(x, y, control);
}

// [[Rcpp::export]]
Rcpp::List
golemDense(const arma::mat& x,
           const arma::mat& y,
           const Rcpp::List control)
{
  return golemCpp(x, y, control);
}
