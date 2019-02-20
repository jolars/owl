#ifndef GOLEM_SOLVERS_
#define GOLEM_SOLVERS_

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

    double L = family->lipschitzConstant(X);
    double t = 1;
    double t_old = t;

    uword i = 0;
    bool accepted = false;

    ConvergenceCheck convergenceCheck{beta, tol};

    while (!accepted && i < max_passes) {
      // gradient
      pseudo_g = family->gradient((X * beta) + intercept, y);
      g = X.t() * pseudo_g;
      g_intercept = mean(pseudo_g);

      // store old values
      intercept_tilde_old = intercept_tilde;
      beta_tilde_old = beta_tilde;

      // update beta and intercept
      beta_tilde = beta - (1.0/L)*g;
      if (fit_intercept)
        intercept_tilde = intercept - (1.0/L)*g_intercept;

      // apply regularization
      beta_tilde = penalty->eval(beta_tilde, L);

      // FISTA step
      t_old = t;
      t = 0.5*(1.0 + std::sqrt(1.0 + 4.0*t_old*t_old));
      beta = beta_tilde + (t_old - 1.0)/t * (beta_tilde - beta_tilde_old);

      if (fit_intercept)
        intercept = intercept_tilde + (t_old - 1.0)/t * (intercept_tilde - intercept_tilde_old);

      accepted = convergenceCheck(beta);

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

#endif /* GOLEM_SOLVERS_ */

