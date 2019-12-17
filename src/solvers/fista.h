#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "../families.h"
#include "../penalties.h"
#include "solver.h"

using namespace arma;
using namespace Rcpp;

class FISTA : public Solver {
private:
  const bool standardize;
  const bool is_sparse;
  const bool diagnostics;
  const uword max_passes;
  const double tol_rel_gap;
  const double tol_infeas;
  const bool line_search;
  const uword line_search_frequency;
  const uword verbosity;

public:
  FISTA(const bool standardize,
        const bool is_sparse,
        const bool diagnostics,
        const uword max_passes,
        const double tol_rel_gap,
        const double tol_infeas,
        const bool line_search,
        const uword line_search_frequency,
        const uword verbosity)
        : standardize(standardize),
          is_sparse(is_sparse),
          diagnostics(diagnostics),
          max_passes(max_passes),
          tol_rel_gap(tol_rel_gap),
          tol_infeas(tol_infeas),
          line_search(line_search),
          line_search_frequency(line_search_frequency),
          verbosity(verbosity) {}

  virtual
  Results
  fit(const mat& x,
      const mat& y,
      const std::unique_ptr<Family>& family,
      const std::unique_ptr<Penalty>& penalty,
      const rowvec& intercept_init,
      const mat& beta_init,
      const bool fit_intercept,
      const double lipschitz_constant,
      const vec& lambda,
      const vec& x_center,
      const vec& x_scale)
  {
    return fitImpl(x,
                   y,
                   family,
                   penalty,
                   intercept_init,
                   beta_init,
                   fit_intercept,
                   lipschitz_constant,
                   lambda,
                   x_center,
                   x_scale);
  }

  virtual
  Results
  fit(const sp_mat& x,
      const mat& y,
      const std::unique_ptr<Family>& family,
      const std::unique_ptr<Penalty>& penalty,
      const rowvec& intercept_init,
      const mat& beta_init,
      const bool fit_intercept,
      const double lipschitz_constant,
      const vec& lambda,
      const vec& x_center,
      const vec& x_scale)
  {
    return fitImpl(x,
                   y,
                   family,
                   penalty,
                   intercept_init,
                   beta_init,
                   fit_intercept,
                   lipschitz_constant,
                   lambda,
                   x_center,
                   x_scale);
  }

  template <typename T>
  Results
  fitImpl(const T& x,
          const mat& y,
          const std::unique_ptr<Family>& family,
          const std::unique_ptr<Penalty>& penalty,
          const rowvec& intercept_init,
          const mat& beta_init,
          const bool fit_intercept,
          const double lipschitz_constant,
          const vec& lambda,
          const vec& x_center,
          const vec& x_scale)
  {
    using namespace arma;

    uword n = y.n_rows;
    uword p = x.n_cols;
    uword m = beta_init.n_cols;

    rowvec intercept(intercept_init);
    rowvec intercept_tilde(intercept_init);
    rowvec intercept_tilde_old(intercept_init);

    mat beta(beta_init);
    mat beta_tilde(beta_init);
    mat beta_tilde_old(beta_init);

    mat lin_pred(n, m);

    mat gradient(p, m, fill::zeros);
    mat pseudo_gradient(n, m, fill::zeros);
    rowvec gradient_intercept(m, fill::zeros);

    double learning_rate = 1.0/lipschitz_constant;

    // line search parameters
    double eta = 0.5;

    // FISTA parameters
    double t = 1;

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

    lin_pred = linearPredictor(x,
                               beta,
                               intercept,
                               x_center,
                               x_scale,
                               standardize);

    // main loop
    uword passes = 0;
    while (passes < max_passes) {
      if (verbosity >= 3)
        Rcout << "pass: " << passes << std::endl;

      ++passes;

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

      if (verbosity >= 3) {
        Rcout << "coefficients:" << std::endl;
        beta.print();
        Rcout << std::endl;

        Rcout << "gradient:" << std::endl;
        gradient.print();
        Rcout << std::endl;
      }

      double primal = f + penalty->primal(beta, lambda);
      double dual = family->dual(y, lin_pred);
      double infeasibility = penalty->infeasibility(gradient, lambda);

      if (verbosity >= 2) {
        Rcout << "primal: " << primal << "\t"
              << "dual: "   << dual << "\t"
              << "infeas: " << infeasibility << std::endl;
      }

      double small = std::sqrt(datum::eps);

      bool accepted =
        (std::abs(primal - dual)/std::max(small, primal) < tol_rel_gap)
        && (infeasibility <= std::max(small, lambda(0)*tol_infeas));

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
      if (line_search && passes % line_search_frequency == 0) {
        while (true) {
          current_line_searches++;

          // Update beta and intercept
          beta_tilde = penalty->eval(beta - learning_rate*gradient,
                                     lambda,
                                     learning_rate);

          vec d = vectorise(beta_tilde - beta);

          if (fit_intercept)
            intercept_tilde = intercept - learning_rate*gradient_intercept;

          lin_pred = linearPredictor(x,
                                     beta_tilde,
                                     intercept_tilde,
                                     x_center,
                                     x_scale,
                                     standardize);

          f = family->primal(y, lin_pred);
          double q = f_old
            + dot(d, vectorise(gradient))
            + (1.0/(2*learning_rate))*accu(square(d));

          if (q >= f*(1 - 1e-12)) {
            break;
          } else {
            learning_rate *= eta;
          }

          Rcpp::checkUserInterrupt();
        }
      } else {
        // Update beta and intercept
        beta_tilde = penalty->eval(beta - learning_rate*gradient,
                                   lambda,
                                   learning_rate);

        if (fit_intercept)
          intercept_tilde = intercept - learning_rate*gradient_intercept;

        lin_pred = linearPredictor(x,
                                   beta_tilde,
                                   intercept_tilde,
                                   x_center,
                                   x_scale,
                                   standardize);
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
