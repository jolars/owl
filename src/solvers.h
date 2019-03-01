#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "penalties.h"
#include "families.h"
#include "utils.h"

struct Results {
  Results(arma::rowvec intercept,
          arma::mat beta,
          arma::uword passes,
          std::vector<double> loss,
          std::vector<double> time)
          : intercept(intercept),
            beta(beta),
            passes(passes),
            loss(loss),
            time(time) {}

  arma::rowvec intercept;
  arma::mat beta;
  arma::uword passes;
  std::vector<double> loss;
  std::vector<double> time;
};

class Solver {
protected:
  bool diagnostics;
  std::vector<double> loss;
  std::vector<double> time;
  arma::uword path_iter = 0;

  arma::rowvec intercept;
  arma::mat beta;

public:
  virtual
  Results
  fit(const arma::mat& x,
      const arma::mat& y,
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
  FISTA(arma::rowvec&& intercept_init,
        arma::mat&& beta_init,
        const Rcpp::List& args)
  {
    using Rcpp::as;

    intercept = std::move(intercept_init);
    beta = std::move(beta_init);

    max_passes  = as<arma::uword>(args["max_passes"]);
    tol         = as<double>(args["tol"]);
    diagnostics = as<bool>(args["diagnostics"]);
  }

  Results
  fit(const arma::mat& x,
      const arma::mat& y,
      const std::unique_ptr<Family>& family,
      const std::unique_ptr<Penalty>& penalty,
      const bool fit_intercept)
  {
    using namespace arma;

    uword n = x.n_rows;
    uword p = x.n_cols;
    uword n_responses = y.n_cols;

    rowvec intercept_tilde(intercept);
    rowvec intercept_tilde_old(intercept_tilde);

    mat beta_tilde(beta);
    mat beta_tilde_old(beta_tilde);

    mat lin_pred(n, n_responses);

    mat g(p, n_responses, fill::zeros);
    mat pseudo_g(g);
    rowvec g_intercept(n_responses, fill::zeros);

    double L = family->lipschitzConstant(x, fit_intercept);
    //double L = 1.0;
    double t = 1;

    uword i = 0;
    bool accepted = false;

    ConvergenceCheck convergenceCheck{intercept, beta, tol};

    // diagnostics
    wall_clock timer;

    if (diagnostics) {
      loss.reserve(max_passes);
      time.reserve(max_passes);
      timer.tic();
    }

    lin_pred = x*beta;
    if (fit_intercept)
      lin_pred.each_row() += intercept;

    while (!accepted && i < max_passes) {
      // loss and gradient
      double f = family->loss(lin_pred, y);
      pseudo_g = family->gradient(lin_pred, y);
      g = x.t() * pseudo_g;

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
        mat d = beta_tilde - beta;
        if (fit_intercept)
          intercept_tilde = intercept - (1.0/L)*g_intercept;

        lin_pred = x*beta_tilde;
        if (fit_intercept)
          lin_pred.each_row() += intercept;

        f = family->loss(lin_pred, y);
        mat q = f_old + d.t()*g + 0.5*L*d.t()*d;

        if (any(q.diag() >= f*(1 - 1e-12)))
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

      lin_pred = x*beta;
      if (fit_intercept)
        lin_pred.each_row() += intercept;

      if (diagnostics) {
        time.push_back(timer.toc());
        loss.push_back(family->loss(lin_pred, y) + penalty->primal(beta));
      }

      accepted = convergenceCheck(intercept, beta);

      i++;
    }

    Results res{intercept, beta, i, loss, time};

    return res;
  }
};

