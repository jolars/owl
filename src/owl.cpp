#include <RcppArmadillo.h>
#include <memory>
#include "solvers.h"
#include "penalties.h"
#include "families.h"
#include "screening_rules.h"
#include "lipschitzConstant.h"

template <typename T>
Rcpp::List
owlCpp(const T& x,
       const arma::vec& y,
       const Rcpp::List control)
{
  using namespace arma;
  using Rcpp::as;
  using Rcpp::Named;
  using Rcpp::wrap;

  // auto n = x.n_rows;
  auto p = x.n_cols;
  auto n = x.n_rows;

  // parameter packs for penalty and solver
  auto penalty_args = as<Rcpp::List>(control["penalty"]);
  auto solver_args = as<Rcpp::List>(control["solver"]);

  auto max_passes = as<uword>(control["max_passes"]);
  auto diagnostics = as<bool>(control["diagnostics"]);

  auto groups = as<Rcpp::List>(control["groups"]);

  auto family_args = as<Rcpp::List>(control["family"]);
  auto fit_intercept = as<bool>(control["fit_intercept"]);
  auto screening_rule = as<std::string>(control["screening_rule"]);
  auto sigma_type = as<std::string>(control["sigma_type"]);

  auto standardize_features = as<bool>(control["standardize_features"]);
  auto is_sparse = as<bool>(control["is_sparse"]);
  auto n_sigma = as<uword>(control["n_sigma"]);

  auto lambda = as<vec>(control["lambda"]);
  auto sigma  = as<vec>(control["sigma"]);

  // get scaled vector of feature matrix centers for use in sparse fitting
  vec x_center        = as<vec>(control["x_center"]);
  vec x_scale         = as<vec>(control["x_scale"]);
  vec x_scaled_center = x_center/x_scale;

  auto family_choice = as<std::string>(family_args["name"]);

  // setup family and response
  auto family = setupFamily(family_choice,
                            fit_intercept,
                            standardize_features);

  auto penalty = setupPenalty(penalty_args, groups);

  cube betas(p, 1, n_sigma, fill::zeros);
  cube intercepts(1, 1, n_sigma);

  // initialize estimates
  auto intercept = as<double>(control["intercept_init"]);
  auto beta      = as<vec>(control["beta_init"]);

  if (fit_intercept)
    intercept = family->fitNullModel(y);

  double intercept_prev = intercept;
  vec beta_prev(p, fill::zeros);

  uvec passes(n_sigma);
  std::vector<std::vector<double>> primals;
  std::vector<std::vector<double>> duals;
  std::vector<std::vector<double>> timings;
  std::vector<std::vector<double>> infeasibilities;
  std::vector<std::vector<unsigned>> line_searches;
  uvec violations(n_sigma, fill::zeros);

  vec linear_predictor_prev(n);
  vec gradient_prev(p);
  vec pseudo_gradient_prev(n);

  // sets of active predictors
  field<uvec> active_sets(n_sigma);
  uvec active_set;

  auto solver = setupSolver(as<std::string>(solver_args["name"]),
                            standardize_features,
                            is_sparse,
                            diagnostics,
                            max_passes,
                            solver_args);

  vec x_colnorms(p);

  if (screening_rule == "safe") {
    if (standardize_features && is_sparse) {
      for (uword j = 0; j < p; ++j)
        x_colnorms(j) = norm(x.col(j) - x_center(j), 2);

    } else {
      for (uword j = 0; j < p; ++j)
        x_colnorms(j) = norm(x.col(j), 2);
    }
  }

  std::vector<double> lipschitz_constants;
  double lipschitz_constant{1};

  if (screening_rule == "none") {
    lipschitz_constant = lipschitzConstant(x,
                                           x_center,
                                           x_scale,
                                           standardize_features,
                                           family_choice);
  }

  Results res;

  for (uword k = 0; k < n_sigma; ++k) {

    if (screening_rule == "none") {

      active_set = regspace<uvec>(0, p-1);

    } else {

      if (sigma_type == "sequence" && k == 0) {
        // no predictors active at first step (except intercept)

        active_set.set_size(0);
        beta.zeros();

      } else {

        // NOTE(JL): the screening rules should probably not be used if
        // the coefficients from the previous fit are already very dense

        linear_predictor_prev = linearPredictor(x,
                                                beta_prev,
                                                intercept_prev,
                                                x_center,
                                                x_scale,
                                                standardize_features);

        pseudo_gradient_prev = family->pseudoGradient(y, linear_predictor_prev);
        gradient_prev = x.t() * pseudo_gradient_prev;

        active_set = penalty->activeSet(family,
                                        y,
                                        gradient_prev,
                                        pseudo_gradient_prev,
                                        x_colnorms,
                                        lambda*sigma(k),
                                        lambda*sigma(k-1),
                                        screening_rule);
      }
    }

    if (active_set.n_elem == 0) {
      // null (intercept only) model

      if (fit_intercept)
        intercept = family->fitNullModel(y);

      beta.zeros();

      passes(k) = 0;

      if (diagnostics) {
        primals.emplace_back(0);
        duals.emplace_back(0);
        infeasibilities.emplace_back(0);
        timings.emplace_back(0);
        line_searches.emplace_back(0);
      }

    } else if (active_set.n_elem == p) {
      // all features active

      res = solver->fit(x,
                        y,
                        family,
                        penalty,
                        intercept,
                        beta,
                        fit_intercept,
                        lipschitz_constant,
                        lambda*sigma(k),
                        x_center,
                        x_scale);

      passes(k) = res.passes;
      beta = res.beta;
      intercept = res.intercept;

      if (diagnostics) {
        primals.push_back(res.primals);
        duals.push_back(res.duals);
        infeasibilities.push_back(res.infeasibilities);
        timings.push_back(res.time);
        line_searches.emplace_back(res.line_searches);
        lipschitz_constants.emplace_back(lipschitz_constant);
      }

    } else {

      bool kkt_violation = true;

      do {
        T x_subset = matrixSubset(x, active_set);

        // compute Lipschitz estimate on x subset
        // NOTE(JL): is it worth it to do this each time?
        lipschitz_constant = lipschitzConstant(x_subset,
                                               x_center(active_set),
                                               x_scale(active_set),
                                               standardize_features,
                                               family_choice);

        res = solver->fit(x_subset,
                          y,
                          family,
                          penalty,
                          intercept,
                          beta(active_set),
                          fit_intercept,
                          lipschitz_constant,
                          lambda.head(active_set.n_elem)*sigma(k),
                          x_center(active_set),
                          x_scale(active_set));

        beta(active_set) = res.beta;
        intercept = res.intercept;
        passes(k) = res.passes;

        linear_predictor_prev = linearPredictor(x,
                                                beta,
                                                intercept,
                                                x_center,
                                                x_scale,
                                                standardize_features);

        pseudo_gradient_prev = family->pseudoGradient(y, linear_predictor_prev);
        gradient_prev = x.t() * pseudo_gradient_prev;

        uvec check_failures =
          penalty->kktCheck(gradient_prev, beta, lambda*sigma(k), 1e-6);

        kkt_violation = check_failures.n_elem > 0;
        violations(k) += check_failures.n_elem;

        active_set = setUnion(check_failures, active_set);

        Rcpp::checkUserInterrupt();

      } while (kkt_violation);

      if (diagnostics) {
        primals.push_back(res.primals);
        duals.push_back(res.duals);
        infeasibilities.push_back(res.infeasibilities);
        timings.push_back(res.time);
        line_searches.push_back(res.line_searches);
        lipschitz_constants.emplace_back(lipschitz_constant);
      }
    }

    // store coefficients and intercept
    betas.slice(k) = beta;
    intercepts.slice(k) = intercept;
    intercept_prev = intercept;
    beta_prev = beta;

    active_sets(k) = active_set;

    Rcpp::checkUserInterrupt();
  }

  return Rcpp::List::create(
    Named("intercept")           = wrap(intercepts),
    Named("beta")                = wrap(betas),
    Named("active_sets")         = wrap(active_sets),
    Named("passes")              = passes,
    Named("primals")             = wrap(primals),
    Named("duals")               = wrap(duals),
    Named("infeasibilities")     = wrap(infeasibilities),
    Named("time")                = wrap(timings),
    Named("line_searches")       = wrap(line_searches),
    Named("violations")          = wrap(violations),
    Named("lipschitz_constants") = wrap(lipschitz_constants)
  );
}


// [[Rcpp::export]]
Rcpp::List
owlSparse(const arma::sp_mat& x,
          const arma::vec& y,
          const Rcpp::List control)
{
  return owlCpp(x, y, control);
}

// [[Rcpp::export]]
Rcpp::List
owlDense(const arma::mat& x,
         const arma::vec& y,
         const Rcpp::List control)
{
  return owlCpp(x, y, control);
}
