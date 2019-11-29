#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "penalties.h"
#include "families.h"
#include "utils.h"

struct Results {
  Results(double intercept,
          arma::vec beta,
          arma::vec gradient,
          arma::uword passes,
          std::vector<double> primals,
          std::vector<double> duals,
          std::vector<double> infeasibilities,
          std::vector<double> time)
          : intercept(intercept),
            beta(beta),
            gradient(gradient),
            passes(passes),
            primals(primals),
            duals(duals),
            infeasibilities(infeasibilities),
            time(time) {}

  double intercept;
  arma::vec beta;
  arma::vec gradient;
  arma::uword passes;
  std::vector<double> primals;
  std::vector<double> duals;
  std::vector<double> infeasibilities;
  std::vector<double> time;
};

class Solver {
protected:
  bool diagnostics;
  double intercept;
  arma::vec beta;
};

class FISTA : public Solver {
private:
  const bool standardize;
  arma::vec x_scaled_center;
  const bool is_sparse;
  arma::uword max_passes;
  const double eta = 2.0;
  double tol_rel_gap = 1e-6;
  double tol_infeas = 1e-6;

public:
  FISTA(const double intercept_init,
        const arma::vec& beta_init,
        const bool standardize,
        const arma::vec x_scaled_center,
        const bool is_sparse,
        const Rcpp::List& args)
        : standardize(standardize),
          x_scaled_center(x_scaled_center),
          is_sparse(is_sparse)
  {
    using Rcpp::as;

    intercept   = intercept_init;
    beta        = beta_init;

    max_passes  = as<arma::uword>(args["max_passes"]);
    diagnostics = as<bool>(args["diagnostics"]);
    tol_rel_gap = as<double>(args["tol_rel_gap"]);
    tol_infeas  = as<double>(args["tol_infeas"]);
  }

  template <typename T>
  Results
  fit(const T& x,
      const arma::vec& y,
      const std::unique_ptr<Family>& family,
      const std::unique_ptr<Penalty>& penalty,
      const bool fit_intercept,
      double L,
      const arma::vec& lambda)
  {
    using namespace arma;

    uword n = y.n_elem;
    uword p = x.n_cols;

    double intercept_tilde = intercept;
    double intercept_tilde_old = intercept_tilde;

    vec beta_tilde(beta);
    vec beta_tilde_old(beta_tilde);

    vec lin_pred(n);

    vec gradient(p, fill::zeros);
    vec pseudo_gradient(gradient);
    double gradient_intercept = 0.0;

    double t = 1;

    uword passes = 0;
    bool accepted = false;

    double mod_tol_infeas{tol_infeas};

    mod_tol_infeas *= lambda(0);

    if (mod_tol_infeas < std::sqrt(datum::eps))
      mod_tol_infeas = std::sqrt(datum::eps);

    // diagnostics
    wall_clock timer;
    std::vector<double> primals;
    std::vector<double> duals;
    std::vector<double> infeasibilities;
    std::vector<double> time;

    if (diagnostics) {
      primals.reserve(max_passes);
      duals.reserve(max_passes);
      infeasibilities.reserve(max_passes);
      time.reserve(max_passes);
      timer.tic();
    }

    if (standardize && is_sparse)
      lin_pred = x*beta - arma::dot(x_scaled_center, beta);
    else
      lin_pred = x*beta;

    if (fit_intercept)
      lin_pred += intercept;

    // main loop
    while (!accepted && passes < max_passes) {
      ++passes;
      // gradient
      double f = family->primal(y, lin_pred);
      pseudo_gradient = family->pseudoGradient(y, lin_pred);
      gradient = x.t()*pseudo_gradient;

      // adjust gradient if sparse and standardizing
      if (is_sparse && standardize) {
        for (uword j = 0; j < p; ++j)
          gradient(j) -= accu(x_scaled_center(j)*pseudo_gradient);
      }

      if (fit_intercept)
        gradient_intercept = mean(pseudo_gradient);

      double primal = f + penalty->primal(beta, lambda);
      double dual = family->dual(y, lin_pred);
      double infeasibility = penalty->infeasibility(gradient, lambda);

      accepted = (std::abs(primal - dual)/std::max(1.0, primal) < tol_rel_gap)
                  && (infeasibility <= mod_tol_infeas);

      if (diagnostics) {
        time.push_back(timer.toc());
        primals.push_back(primal);
        infeasibilities.push_back(infeasibility);
        duals.push_back(dual);
      }

      if (accepted)
        break;

      beta_tilde_old = beta_tilde;

      if (fit_intercept)
        intercept_tilde_old = intercept_tilde;

      double f_old = f;
      double t_old = t;

      // // Lipschitz search
      while (true) {
        // Update beta and intercept
        beta_tilde = penalty->eval(beta - (1.0/L)*gradient, lambda, 1.0/L);

        vec d = beta_tilde - beta;

        if (fit_intercept)
          intercept_tilde = intercept - (1.0/L)*gradient_intercept;

        if (standardize && is_sparse)
          lin_pred = x*beta_tilde - arma::dot(x_scaled_center, beta_tilde);
        else
          lin_pred = x*beta_tilde;

        if (fit_intercept)
          lin_pred += intercept_tilde;

        f = family->primal(y, lin_pred);
        mat q = f_old + d.t()*gradient + 0.5*L*d.t()*d;

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

      if (standardize && is_sparse)
        lin_pred = x*beta - arma::dot(x_scaled_center, beta);
      else
        lin_pred = x*beta;

      if (fit_intercept)
        lin_pred += intercept;

      if (passes % 100 == 0)
        Rcpp::checkUserInterrupt();
    }

    Results res{intercept,
                beta,
                gradient,
                passes,
                primals,
                duals,
                infeasibilities,
                time};

    return res;
  }
};

