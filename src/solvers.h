#ifndef GOLEM_SOLVERS_
#define GOLEM_SOLVERS_

#include "penalties.h"
#include "families.h"
#include "utils.h"

class Solver {
public:
  virtual Rcpp::List fit() = 0;
};

class FISTA : public Solver {
public:
  FISTA(const arma::mat& X,
        const arma::vec& y,
        const arma::vec& lambda,
        const std::unique_ptr<Family>& family,
        const std::unique_ptr<Penalty>& penalty)
        : X(X), y(y), lambda(lambda), family(family), penalty(penalty) {}

  Rcpp::List
  fit()
  {
    using namespace arma;

    uword n = X.n_rows;
    uword p = X.n_cols;

    vec beta(p, fill::zeros);
    vec beta_tilde(beta);
    vec beta_tilde_old(beta_tilde);
    vec g(p, fill::zeros);

    double L = family->lipschitzConstant();
    double t = 1.0;
    double t_old = t;

    uword i = 0;
    bool accepted = false;

    ConvergenceCheck convergenceCheck{beta, tol};

    while (!accepted && i < max_iter) {
      // gradient
      g = family->gradient(X * beta);

      beta_tilde_old = beta_tilde;
      beta_tilde = penalty->eval(beta - (1.0/L)*g, L);

      t_old = t;
      t = 0.5*(1.0 + std::sqrt(1.0 + 4.0*t_old*t_old));
      beta = beta_tilde + (t_old - 1.0)/t * (beta_tilde - beta_tilde_old);

      accepted = convergenceCheck(beta);

      i++;
    }

    return Rcpp::List::create(
      Rcpp::Named("x")          = Rcpp::wrap(beta),
      Rcpp::Named("iter")       = i,
      Rcpp::Named("lipschitz")  = L
    );
  }

private:
  const arma::mat& X;
  const arma::vec& y;
  const arma::vec& lambda;
  const std::unique_ptr<Family>& family;
  const std::unique_ptr<Penalty>& penalty;
  const arma::uword max_iter = 1e4;
  double tol = 1e-4;
};

#endif /* GOLEM_SOLVERS_ */

