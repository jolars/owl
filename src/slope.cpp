#include <RcppArmadillo.h>
#include "proxes.h"
#include "families.h"

std::unique_ptr<Prox>
setupProx(const std::string& prox_choice,
          const arma::vec& lambda)
{
  return std::unique_ptr<SLOPE>(new SLOPE{lambda});
}

std::unique_ptr<Family>
setupFamily(const std::string& family_choice,
            const arma::mat& X,
            const arma::vec& y)
{
  return std::unique_ptr<Gaussian>(new Gaussian(X, y));
}

// [[Rcpp::export]]
Rcpp::List
slope_solver(const arma::mat X,
             const arma::vec y,
             const arma::vec lambda,
             arma::uword max_iter = 1e4,
             arma::uword grad_iter = 20,
             arma::uword opt_iter = 1,
             double tol_infeas = 1e-6,
             double tol_rel_gap = 1e-6)
{
  using namespace arma;

  uword n = X.n_cols;

  // select prox
  auto prox = setupProx("slope", lambda);

  // select family
  auto family = setupFamily("gaussian", X, y);

  // lower bound on Lipschitz constant
  vec x = Rcpp::rnorm(n);
  x = x / norm(x, 2);
  x = X.t() * X * x;
  double L = norm(x, 2);

  // initialize parameters and iterates
  //vec x(n, fill::zeros);
  double t = 1.0;
  double t_prev = 1.0;

  double eta = 2.0;

  vec beta(n, fill::zeros);
  vec beta_tilde(n, fill::zeros);

  mat X_beta_tilde = X*beta_tilde;
  mat X_beta_tilde_prev(X_beta_tilde);

  double f_prev = datum::inf;
  uword iter = 0;
  bool optimal = false;
  double infeas = 0.0;
  double obj_primal = 0.0;
  double obj_dual = 0.0;

  while (iter <= max_iter) {
    vec lin_pred(n);

    // compute gradient at f(beta)
    if (iter & grad_iter == 0)
      lin_pred = X*beta;
    else
      lin_pred = (X_beta_tilde + ((t_prev - 1.0)/t) * (X_beta_tilde - X_beta_tilde_prev));

    vec g = family->gradient(lin_pred);
    double f = family->loss(lin_pred);

    iter++;

    if (iter % opt_iter == 0) {
      // compute dual, check infeasibility and gap
      vec g_sorted = sort(abs(g), "descending");
      vec beta_sorted = sort(abs(beta), "descending");
      infeas = std::max(cumsum(g_sorted - lambda).max(), 0.0);

      // compute primal and dual objective
      obj_primal =  f + dot(lambda, beta_sorted);
      obj_dual   = -f - dot(lin_pred - y, y);

      // check primal-dual gap
      if ((std::abs(obj_primal - obj_dual)/std::max(1.0, obj_primal) < tol_rel_gap) &&
        (infeas < tol_infeas*lambda(0))) {
          optimal = true;
      }
    }

    // stopping criteria
    if (optimal)
      break;

    // store previous values
    X_beta_tilde_prev = X_beta_tilde;
    vec beta_tilde_prev = beta_tilde;
    f_prev = f;
    t_prev = t;

    // Lipschitz search
    while (true) {
      // compute prox mapping
      beta_tilde = prox->eval(beta - (1.0/L)*g, L);

      vec d = beta_tilde - beta;
      X_beta_tilde = X*beta_tilde;
      f = family->loss(X_beta_tilde);
      double q = f_prev + dot(d, g) + (L/2.0)*dot(d, d);

      if (q < f*(1.0 - 1e-12))
        L *= eta;
      else
        break;
    }

    // update
    t = (1.0 + std::sqrt(1.0 + 4.0*t*t))/2.0;
    beta = beta_tilde + ((t_prev - 1.0)/t) * (beta_tilde - beta_tilde_prev);

  }

  return Rcpp::List::create(
    Rcpp::Named("x")          = Rcpp::wrap(beta),
    Rcpp::Named("optimal")    = optimal,
    Rcpp::Named("iter")       = iter,
    Rcpp::Named("infeas")     = infeas,
    Rcpp::Named("obj_primal") = obj_primal,
    Rcpp::Named("obj_dual")   = obj_dual,
    Rcpp::Named("lipschitz")  = L
  );
}



