#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "penalties.h"
#include "families.h"
#include "utils.h"

struct Results {
  double intercept;
  arma::vec beta;
  arma::uword passes;
  std::vector<double> primals;
  std::vector<double> duals;
  std::vector<double> infeasibilities;
  std::vector<double> time;
  std::vector<unsigned> line_searches;

  Results() {}

  Results(double intercept,
          arma::vec beta,
          arma::uword passes,
          std::vector<double> primals,
          std::vector<double> duals,
          std::vector<double> infeasibilities,
          std::vector<double> time,
          std::vector<unsigned> line_searches)
          : intercept(intercept),
            beta(beta),
            passes(passes),
            primals(primals),
            duals(duals),
            infeasibilities(infeasibilities),
            time(time),
            line_searches(line_searches) {}
};

class Solver {
protected:
  bool diagnostics;
};

class FISTA : public Solver {
private:
  const bool standardize;
  const bool is_sparse;
  arma::uword max_passes;
  double tol_rel_gap = 1e-6;
  double tol_infeas = 1e-6;

public:
  FISTA(const bool standardize,
        const bool is_sparse,
        const Rcpp::List& args)
        : standardize(standardize),
          is_sparse(is_sparse)
  {
    using Rcpp::as;

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
      const double intercept_init,
      const arma::vec& beta_init,
      const bool fit_intercept,
      double lipschitz_constant,
      const arma::vec& lambda,
      const arma::vec& x_center,
      const arma::vec& x_scale)
  {
    using namespace arma;

    uword n = y.n_elem;
    uword p = x.n_cols;

    double intercept = intercept_init;
    double intercept_tilde = intercept;
    double intercept_tilde_old = intercept_tilde;

    vec beta(beta_init);
    vec beta_tilde(beta);
    vec beta_tilde_old(beta_tilde);

    vec lin_pred(n);

    vec gradient(p, fill::zeros);
    vec pseudo_gradient(n, fill::zeros);
    double gradient_intercept = 0.0;

    double learning_rate = 1.0/lipschitz_constant;
    // double learning_rate = 1;

    // line search parameters
    double eta = 0.5;

    // FISTA parameters
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
    std::vector<unsigned> line_searches;

    if (diagnostics) {
      primals.reserve(max_passes);
      duals.reserve(max_passes);
      infeasibilities.reserve(max_passes);
      time.reserve(max_passes);
      line_searches.reserve(max_passes);
      timer.tic();
    }

    if (standardize && is_sparse)
      lin_pred = x*beta - arma::dot(x_center/x_scale, beta);
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
          gradient(j) -= accu(pseudo_gradient*x_center(j)/x_scale(j));
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

      if (accepted) {
        if (diagnostics)
          line_searches.emplace_back(0);
        break;
      }

      beta_tilde_old = beta_tilde;

      if (fit_intercept)
        intercept_tilde_old = intercept_tilde;

      double f_old = f;
      double t_old = t;

      unsigned current_line_searches = 0;

      // Backtracking line search
      while (true) {
        current_line_searches++;

        // Update beta and intercept
        beta_tilde = penalty->eval(beta - learning_rate*gradient,
                                   lambda,
                                   learning_rate);

        vec d = beta_tilde - beta;

        if (fit_intercept)
          intercept_tilde = intercept - learning_rate*gradient_intercept;

        if (standardize && is_sparse)
          lin_pred = x*beta_tilde - arma::dot(x_center/x_scale, beta_tilde);
        else
          lin_pred = x*beta_tilde;

        if (fit_intercept)
          lin_pred += intercept_tilde;

        f = family->primal(y, lin_pred);
        double q =
          f_old + dot(d, gradient) + (1.0/(2*learning_rate))*accu(square(d));

        if (q >= f*(1 - 1e-12)) {
          break;
        } else {
          learning_rate *= eta;
        }
      }

      if (diagnostics)
        line_searches.emplace_back(current_line_searches);

      // FISTA step
      t = 0.5*(1.0 + std::sqrt(1.0 + 4.0*t_old*t_old));
      beta = beta_tilde + (t_old - 1.0)/t * (beta_tilde - beta_tilde_old);

      if (fit_intercept)
        intercept = intercept_tilde
                    + (t_old - 1.0)/t * (intercept_tilde - intercept_tilde_old);

      lin_pred = linearPredictor(x,
                                 beta,
                                 intercept,
                                 x_center,
                                 x_scale,
                                 standardize);

      if (passes % 100 == 0)
        Rcpp::checkUserInterrupt();
    }

    Results res{intercept,
                beta,
                passes,
                primals,
                duals,
                infeasibilities,
                time,
                line_searches};

    return res;
  }
};

