#include <RcppArmadillo.h>
#include <memory>
#include "solvers.h"
#include "penalties.h"
#include "families.h"
#include "screening_rules.h"
#include "standardize.h"
#include "rescale.h"
#include "regularizationPath.h"

using namespace Rcpp;
using namespace arma;

template <typename T>
List owlCpp(T& x, mat& y, const List control)
{
  auto p = x.n_cols;
  auto n = x.n_rows;
  auto m = y.n_cols;

  auto tol_dev_ratio = as<double>(control["tol_dev_ratio"]);
  auto tol_dev_change = as<double>(control["tol_dev_change"]);
  auto max_variables = as<uword>(control["max_variables"]);

  auto diagnostics = as<bool>(control["diagnostics"]);
  auto verbosity = as<int>(control["verbosity"]);

  // solver arguments
  auto max_passes = as<uword>(control["max_passes"]);
  auto tol_rel_gap = as<double>(control["tol_rel_gap"]);
  auto tol_infeas = as<double>(control["tol_infeas"]);

  auto family_args = as<List>(control["family"]);
  auto fit_intercept = as<bool>(control["fit_intercept"]);
  auto screening = as<bool>(control["screening"]);

  auto standardize_features = as<bool>(control["standardize_features"]);
  auto is_sparse = as<bool>(control["is_sparse"]);

  auto y_center = as<rowvec>(control["y_center"]);
  auto y_scale = as<rowvec>(control["y_scale"]);
  rowvec x_center(p, fill::zeros);
  rowvec x_scale(p, fill::ones);

  if (standardize_features)
    standardize(x, x_center, x_scale);

  auto family_choice = as<std::string>(family_args["name"]);

  auto lambda = as<vec>(control["lambda"]);
  auto sigma  = as<vec>(control["sigma"]);
  auto lambda_type = as<std::string>(control["lambda_type"]);
  auto sigma_type = as<std::string>(control["sigma_type"]);
  auto lambda_min_ratio = as<double>(control["lambda_min_ratio"]);
  auto q = as<double>(control["q"]);
  uword n_sigma = sigma.n_elem;
  double sigma_max = 0;

  regularizationPath(sigma,
                     lambda,
                     sigma_max,
                     x,
                     y,
                     x_center,
                     x_scale,
                     y_scale,
                     standardize_features,
                     lambda_type,
                     sigma_type,
                     lambda_min_ratio,
                     q,
                     family_choice,
                     is_sparse);

  // setup family and penalty
  if (verbosity >= 1)
    Rcpp::Rcout << "setting up family" << std::endl;

  auto family = setupFamily(family_choice);

  if (verbosity >= 1)
    Rcout << "setting up penalty" << std::endl;

  auto penalty = setupPenalty();

  cube betas(p, m, n_sigma, fill::zeros);
  cube intercepts(1, m, n_sigma);

  rowvec intercept(m, fill::zeros);
  mat beta(p, m, fill::zeros);

  uword n_variables = static_cast<uword>(fit_intercept);

  if (verbosity >= 1)
    Rcout << "fitting intercept for null model" << std::endl;

  if (fit_intercept)
    intercept = family->fitNullModel(y, m);

  mat linear_predictor = linearPredictor(x,
                                         beta,
                                         intercept,
                                         x_center,
                                         x_scale,
                                         standardize_features);

  double null_deviance = 2*family->primal(y, linear_predictor);
  std::vector<double> deviances;
  std::vector<double> deviance_ratios;

  rowvec intercept_prev(m, fill::zeros);
  mat beta_prev(p, m, fill::zeros);

  uvec passes(n_sigma);
  std::vector<std::vector<double>> primals;
  std::vector<std::vector<double>> duals;
  std::vector<std::vector<double>> timings;
  std::vector<std::vector<double>> infeasibilities;
  std::vector<std::vector<unsigned>> line_searches;
  uvec violations(n_sigma, fill::zeros);

  mat linear_predictor_prev(n, m);
  mat gradient_prev(p, m);
  mat pseudo_gradient_prev(n, m);

  // sets of active predictors
  field<uvec> active_sets(n_sigma);
  uvec active_set = regspace<uvec>(0, p-1);

  if (verbosity >= 1)
    Rcout << "setting up solver" << std::endl;

  auto solver = setupSolver("fista",
                            standardize_features,
                            is_sparse,
                            diagnostics,
                            max_passes,
                            tol_rel_gap,
                            tol_infeas,
                            verbosity);

  Results res;

  uword k = 0;

  while (k < n_sigma) {

    if (verbosity >= 1)
      Rcout << "penalty: " << k + 1 << std::endl;

    if (screening) {
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

      double sigma_prev = k == 0 ? sigma_max : sigma(k-1);

      active_set = activeSet(y,
                             gradient_prev,
                             lambda*sigma(k),
                             lambda*sigma_prev);
    }

    if (active_set.n_elem == 0) {
      // null (intercept only) model

      if (fit_intercept)
        intercept = family->fitNullModel(y, m);

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
      }

    } else {

      bool kkt_violation = true;

      do {
        if (verbosity > 0) {
          Rcout << "active predictors:" << std::endl;
          active_set.print();
          Rcout << std::endl;
        }

        T x_subset = matrixSubset(x, active_set);

        res = solver->fit(x_subset,
                          y,
                          family,
                          penalty,
                          intercept,
                          beta.rows(active_set),
                          fit_intercept,
                          lambda.head(active_set.n_elem*m)*sigma(k),
                          x_center.cols(active_set),
                          x_scale.cols(active_set));

        beta.rows(active_set) = res.beta;
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

        uvec possible_failures =
          penalty->kktCheck(gradient_prev, beta, lambda*sigma(k), tol_infeas);
        uvec check_failures = setDiff(possible_failures, active_set);

        if (verbosity >= 1) {
          Rcout << "kkt-failures at: " << std::endl;
          check_failures.print();
          Rcout << std::endl;
        }

        kkt_violation = check_failures.n_elem > 0;
        violations(k) += check_failures.n_elem;

        active_set = setUnion(check_failures, active_set);

        checkUserInterrupt();

      } while (kkt_violation);

      if (diagnostics) {
        primals.push_back(res.primals);
        duals.push_back(res.duals);
        infeasibilities.push_back(res.infeasibilities);
        timings.push_back(res.time);
        line_searches.push_back(res.line_searches);
      }
    }

    // store coefficients and intercept
    double deviance = res.deviance;
    double deviance_ratio = 1.0 - deviance/null_deviance;

    deviances.push_back(deviance);
    deviance_ratios.push_back(deviance_ratio);
    betas.slice(k) = beta;
    intercepts.slice(k) = intercept;
    intercept_prev = intercept;
    beta_prev = beta;

    active_sets(k) = active_set;

    if (verbosity >= 1)
      Rcout << "deviance: " << deviance << "\t"
            << "deviance ratio: " << deviance_ratio << std::endl;

    if (k > 0) {
      // stop path if fractional deviance change is small
      double deviance_change =
        std::abs((deviances[k-1] - deviances[k])/deviances[k-1]);

      if (verbosity >= 1)
        Rcout << "deviance change: " << deviance_change << std::endl;

      if (deviance_change < tol_dev_change || deviance_ratio > tol_dev_ratio) {
        k++;
        break;
      }
    }

    n_variables = static_cast<uword>(fit_intercept) + accu(any(beta != 0, 1));

    if (verbosity >= 1)
      Rcout << "number of variables: " << n_variables << std::endl;

    if (n_variables > max_variables) {
      break;
    }

    k++;

    checkUserInterrupt();
  }

  intercepts.resize(1, m, k);
  betas.resize(p, m, k);
  violations.resize(k);
  passes.resize(k);
  sigma.resize(k);
  active_sets = active_sets.rows(0, k-1);

  rescale(intercepts,
          betas,
          x_center,
          x_scale,
          y_center,
          y_scale,
          fit_intercept);

  return List::create(
    Named("intercepts")          = wrap(intercepts),
    Named("betas")               = wrap(betas),
    Named("active_sets")         = wrap(active_sets),
    Named("passes")              = wrap(passes),
    Named("primals")             = wrap(primals),
    Named("duals")               = wrap(duals),
    Named("infeasibilities")     = wrap(infeasibilities),
    Named("time")                = wrap(timings),
    Named("line_searches")       = wrap(line_searches),
    Named("violations")          = wrap(violations),
    Named("deviance_ratio")      = wrap(deviance_ratios),
    Named("null_deviance")       = wrap(null_deviance),
    Named("sigma")               = wrap(sigma),
    Named("lambda")              = wrap(lambda)
  );
}

// [[Rcpp::export]]
Rcpp::List owlSparse(arma::sp_mat x,
                     arma::mat y,
                     const Rcpp::List control)
{
  return owlCpp(x, y, control);
}

// [[Rcpp::export]]
Rcpp::List owlDense(arma::mat x,
                    arma::mat y,
                    const Rcpp::List control)
{
  return owlCpp(x, y, control);
}
