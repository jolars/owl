#include <RcppArmadillo.h>
#include <memory>
#include "results.h"
#include "families/families.h"
#include "screening.h"
#include "standardize.h"
#include "rescale.h"
#include "regularizationPath.h"
#include "kktCheck.h"

using namespace Rcpp;
using namespace arma;

template <typename T>
List owlCpp(T& x, mat& y, const List control)
{
  auto tol_dev_ratio = as<double>(control["tol_dev_ratio"]);
  auto tol_dev_change = as<double>(control["tol_dev_change"]);
  auto max_variables = as<uword>(control["max_variables"]);

  auto diagnostics = as<bool>(control["diagnostics"]);
  auto verbosity = as<uword>(control["verbosity"]);

  // solver arguments
  auto max_passes  = as<uword>(control["max_passes"]);
  auto tol_rel_gap = as<double>(control["tol_rel_gap"]);
  auto tol_infeas  = as<double>(control["tol_infeas"]);
  auto tol_abs     = as<double>(control["tol_abs"]);
  auto tol_rel     = as<double>(control["tol_rel"]);

  auto family_choice = as<std::string>(control["family"]);
  auto intercept = as<bool>(control["fit_intercept"]);
  auto screening = as<bool>(control["screening"]);

  auto n = x.n_rows;
  auto p = x.n_cols;
  auto m = y.n_cols;

  auto center = as<bool>(control["center"]);
  auto scale = as<std::string>(control["scale"]);

  auto y_center = as<rowvec>(control["y_center"]);
  auto y_scale = as<rowvec>(control["y_scale"]);
  rowvec x_center(p, fill::zeros);
  rowvec x_scale(p, fill::ones);

  standardize(x, x_center, x_scale, intercept, center, scale);

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
                     y_scale,
                     lambda_type,
                     sigma_type,
                     lambda_min_ratio,
                     q,
                     family_choice,
                     intercept);

  auto family = setupFamily(family_choice,
                            intercept,
                            diagnostics,
                            max_passes,
                            tol_rel_gap,
                            tol_infeas,
                            tol_abs,
                            tol_rel,
                            verbosity);

  cube betas(p, m, n_sigma, fill::zeros);
  mat beta(p, m, fill::zeros);

  uword n_variables = 0;
  uvec n_unique(n_sigma);

  if (intercept)
    beta.row(0) = family->fitNullModel(y, m);

  mat linear_predictor = x*beta;

  double null_deviance = 2*family->primal(y, linear_predictor);
  std::vector<double> deviances;
  std::vector<double> deviance_ratios;
  double deviance_change{0};

  mat beta_prev(p, m, fill::zeros);

  uvec passes(n_sigma);
  std::vector<std::vector<double>> primals;
  std::vector<std::vector<double>> duals;
  std::vector<std::vector<double>> timings;
  std::vector<unsigned> violations;
  std::vector<std::vector<unsigned>> violation_list;

  mat linear_predictor_prev(n, m);
  mat gradient_prev(p, m);
  mat pseudo_gradient_prev(n, m);

  // sets of active predictors
  field<uvec> active_sets(n_sigma);
  uvec active_set = regspace<uvec>(0, p-1);

  // object for use in ADMM
  double rho = 0.0;
  vec z(p);
  vec u(p);
  vec z_subset(z);
  vec u_subset(u);
  // for gaussian case
  mat xx, xx_subset, L, U, L_subset, U_subset;
  vec xTy, xTy_subset;

  // factorize xx if gaussian
  if (family->name() == "gaussian") {
    // initialize auxiliary variables
    z.zeros();
    u.zeros();

    // precompute x^Ty
    xTy = x.t() * y;

    // precompute X^tX or XX^t (if wide) and factorize
    if (n >= p) {
      xx = x.t() * x;
    } else {
      xx = x*x.t();
    }

    vec eigval = eig_sym(xx);
    rho = std::pow(eigval.max(), 1/3) * std::pow(lambda.max(), 2/3);

    if (n < p)
      xx /= rho;

    xx.diag() += rho;
  }

  bool factorized = false;

  Results res;

  uword k = 0;

  while (k < n_sigma) {

    violations.clear();

    if (screening) {
      // NOTE(JL): the screening rules should probably not be used if
      // the coefficients from the previous fit are already very dense

      gradient_prev = family->gradient(x, y, x*beta_prev);

      double sigma_prev = k == 0 ? sigma_max : sigma(k-1);

      active_set = activeSet(gradient_prev,
                             lambda*sigma(k),
                             lambda*sigma_prev,
                             intercept);
    }

    if (active_set.n_elem == p/m || !screening) {
      // all features active
      // factorize once if fitting all
      if (!factorized && family->name() == "gaussian") {
        L = chol(xx, "lower");
        U = L.t();
        factorized = true;
      }

      res = family->fit(x, y, beta, z, u, L, U, xTy, lambda*sigma(k), rho);
      passes(k) = res.passes;
      beta = res.beta;

      if (diagnostics) {
        primals.push_back(res.primals);
        duals.push_back(res.duals);
        timings.push_back(res.time);
        violation_list.push_back(violations);
      }

    } else {

      bool kkt_violation = true;

      do {
        T x_subset = matrixSubset(x, active_set);

        if (active_set.n_elem == static_cast<uword>(intercept)) {
          // null model
          beta.zeros();

          if (intercept) {
            // intercept-only model
            beta.row(0) = family->fitNullModel(y, m);
          }

          passes(k) = 0;

        } else {

          if (family->name() == "gaussian") {
            if (n >= p) {
              xx_subset = xx(active_set, active_set);
            } else if (n >= active_set.n_elem) {
              xx_subset = x_subset.t()*x_subset;
              xx_subset.diag() += rho;
            } else {
              xx_subset = x_subset*x_subset.t();
              xx_subset /= rho;
              xx_subset.diag() += rho;
            }

            L_subset = chol(xx_subset, "lower");
            U_subset = L_subset.t();

            z_subset = z(active_set);
            u_subset = u(active_set);
            xTy_subset = xTy(active_set);
          }

          uword n_active = (active_set.n_elem - static_cast<uword>(intercept))*m;

          res = family->fit(x_subset,
                            y,
                            beta.rows(active_set),
                            z_subset,
                            u_subset,
                            L_subset,
                            U_subset,
                            xTy_subset,
                            lambda.head(n_active)*sigma(k),
                            rho);

          beta.rows(active_set) = res.beta;
          passes(k) = res.passes;
        }

        gradient_prev = family->gradient(x, y, x*beta);

        uvec possible_failures =
          kktCheck(gradient_prev, beta, lambda*sigma(k), tol_infeas, intercept);
        uvec check_failures = setDiff(possible_failures, active_set);

        if (verbosity >= 2) {
          Rcout << "kkt-failures at: " << std::endl;
          check_failures.print();
          Rcout << std::endl;
        }

        kkt_violation = check_failures.n_elem > 0;

        if (diagnostics)
          violations.push_back(check_failures.n_elem);

        active_set = setUnion(check_failures, active_set);

        checkUserInterrupt();

      } while (kkt_violation);

      if (diagnostics) {
        primals.push_back(res.primals);
        duals.push_back(res.duals);
        timings.push_back(res.time);
        violation_list.push_back(violations);
      }
    }

    // store coefficients and intercept
    double deviance = res.deviance;
    double deviance_ratio = 1.0 - deviance/null_deviance;
    deviances.push_back(deviance);
    deviance_ratios.push_back(deviance_ratio);

    if (k > 0) {
      deviance_change =
        std::abs((deviances[k-1] - deviance)/deviances[k-1]);
    }

    betas.slice(k) = beta;
    beta_prev = beta;

    active_sets(k) = active_set;
    uword n_coefs = accu(any(beta != 0, 1));
    n_variables = n_coefs;
    n_unique(k) = unique(abs(nonzeros(beta))).eval().n_elem;

    if (verbosity >= 1)
      Rcout << "penalty: "      << k
            << ", dev: "        << deviance
            << ", dev ratio: "  << deviance_ratio
            << ", dev change: " << deviance_change
            << ", n var: "      << n_variables
            << ", n unique: "   << n_unique(k)
            << std::endl;

    if (n_coefs > 0 && k > 0) {
      // stop path if fractional deviance change is small
      if (deviance_change < tol_dev_change || deviance_ratio > tol_dev_ratio) {
        k++;
        break;
      }
    }

    if (n_unique(k) > max_variables)
      break;

    k++;

    checkUserInterrupt();
  }

  betas.resize(p, m, k);
  passes.resize(k);
  sigma.resize(k);
  n_unique.resize(k);
  active_sets = active_sets.rows(0, std::max(static_cast<int>(k-1), 0));

  rescale(betas,
          x_center,
          x_scale,
          y_center,
          y_scale,
          intercept);

  // standardize lambda
  lambda /= n;

  return List::create(
    Named("betas")               = wrap(betas),
    Named("active_sets")         = wrap(active_sets),
    Named("passes")              = wrap(passes),
    Named("primals")             = wrap(primals),
    Named("duals")               = wrap(duals),
    Named("time")                = wrap(timings),
    Named("n_unique")            = wrap(n_unique),
    Named("violations")          = wrap(violation_list),
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

