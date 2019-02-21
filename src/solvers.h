#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "penalties.h"
#include "families.h"
#include "utils.h"

class Solver {
public:
  virtual
  Rcpp::List
  fit(const arma::mat& X,
      const arma::vec& y,
      const std::unique_ptr<Family>& family,
      const std::unique_ptr<Penalty>& penalty,
      const bool fit_intercept) = 0;
};

class FISTA : public Solver {
private:
  arma::uword max_passes;
  double tol;
  const double eta = 2.0;

public:
  FISTA(const Rcpp::List& args)
  {
    using Rcpp::as;

    max_passes = as<arma::uword>(args["max_passes"]);
    tol = as<double>(args["tol"]);
  }

  Rcpp::List
  fit(const arma::mat& X,
      const arma::vec& y,
      const std::unique_ptr<Family>& family,
      const std::unique_ptr<Penalty>& penalty,
      const bool fit_intercept)
  {
    using namespace arma;

    uword n = X.n_rows;
    uword p = X.n_cols;

    double intercept = 0;
    double intercept_tilde = 0;
    double intercept_tilde_old = 0;

    vec beta(p, fill::zeros);
    vec beta_tilde(beta);
    vec beta_tilde_old(beta_tilde);

    vec g(p, fill::zeros);
    vec pseudo_g{g};
    double g_intercept = 0;

    //double L = family->lipschitzConstant(X);
    double L = 1.0;
    double t = 1;

    uword i = 0;
    bool accepted = false;

    ConvergenceCheck convergenceCheck{intercept, beta, tol};

    while (!accepted && i < max_passes) {
      // loss and gradient
      double f = family->loss(X*beta + intercept, y);
      pseudo_g = family->gradient(X*beta + intercept, y);
      g = X.t() * pseudo_g;

      if (fit_intercept) {
        g_intercept = mean(pseudo_g);
        intercept_tilde_old = intercept_tilde;
      }

      beta_tilde_old = beta_tilde;
      double f_old = f;
      double t_old = t;
      beta_tilde_old = beta_tilde;

      // Lipschitz search
      while (true) {
        // Update beta and intercept
        beta_tilde = penalty->eval(beta - (1.0/L)*g, L);
        vec d = beta_tilde - beta;
        if (fit_intercept)
          intercept_tilde = intercept - (1.0/L)*g_intercept;

        f = family->loss(X*beta_tilde + intercept_tilde, y);
        double q = f_old + dot(d, g) + 0.5*L*dot(d, d);

        if (q >= f*(1 - 1e-12))
          break;
        else
          L *= eta;
      }

      // FISTA step
      t = 0.5*(1.0 + std::sqrt(1.0 + 4.0*t_old*t_old));
      beta = beta_tilde + (t_old - 1.0)/t * (beta_tilde - beta_tilde_old);

      if (fit_intercept)
        intercept = intercept_tilde
                    + (t_old - 1.0)/t * (intercept_tilde - intercept_tilde_old);

      accepted = convergenceCheck(intercept, beta);

      i++;
    }

    return Rcpp::List::create(
      Rcpp::Named("intercept")  = intercept,
      Rcpp::Named("beta")       = Rcpp::wrap(beta),
      Rcpp::Named("passes")     = i,
      Rcpp::Named("lipschitz")  = L
    );
  }
};

