#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "families.h"
#include "utils.h"
#include "prox.h"
#include "infeasibility.h"

using namespace arma;
using namespace Rcpp;

struct Results {
  rowvec intercept;
  mat beta;
  uword passes;
  std::vector<double> primals;
  std::vector<double> duals;
  std::vector<double> infeasibilities;
  std::vector<double> time;
  std::vector<unsigned> line_searches;
  double deviance;

  Results() {}

  Results(rowvec intercept,
          mat beta,
          uword passes,
          std::vector<double> primals,
          std::vector<double> duals,
          std::vector<double> infeasibilities,
          std::vector<double> time,
          std::vector<unsigned> line_searches,
          double deviance)
          : intercept(intercept),
            beta(beta),
            passes(passes),
            primals(primals),
            duals(duals),
            infeasibilities(infeasibilities),
            time(time),
            line_searches(line_searches),
            deviance(deviance) {}
};

class Solver {
private:
  const bool standardize;
  const bool is_sparse;
  const bool diagnostics;
  const uword max_passes;
  const double tol_rel_gap;
  const double tol_infeas;
  const uword verbosity;

public:
  Solver(const bool standardize,
         const bool is_sparse,
         const bool diagnostics,
         const uword max_passes,
         const double tol_rel_gap,
         const double tol_infeas,
         const uword verbosity)
    : standardize(standardize),
      is_sparse(is_sparse),
      diagnostics(diagnostics),
      max_passes(max_passes),
      tol_rel_gap(tol_rel_gap),
      tol_infeas(tol_infeas),
      verbosity(verbosity) {}

  template <typename T>
  Results fit(const T& x,
              const mat& y,
              const std::unique_ptr<Family>& family,
              const rowvec& intercept_init,
              const mat& beta_init,
              const bool fit_intercept,
              vec lambda,
              rowvec x_center,
              rowvec x_scale)
  {
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

    double learning_rate = 1.0;

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

    // main loop
    uword passes = 0;
    while (passes < max_passes) {
      if (verbosity >= 3)
        Rcout << "pass: " << passes + 1 << std::endl;

      ++passes;

      linearPredictor(lin_pred,
                      x,
                      beta,
                      intercept,
                      x_center,
                      x_scale,
                      fit_intercept,
                      standardize);

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

      double primal = f + dot(sort(abs(vectorise(beta)), "descending"), lambda);
      double dual = family->dual(y, lin_pred);
      double infeasibility = Infeasibility(gradient, lambda);

      if (verbosity >= 3) {
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
      while (true) {
        current_line_searches++;

        // Update beta and intercept
        beta_tilde = prox(beta - learning_rate*gradient, lambda*learning_rate);

        if (fit_intercept)
          intercept_tilde = intercept - learning_rate*gradient_intercept;

        vec d = vectorise(beta_tilde - beta);

        linearPredictor(lin_pred,
                        x,
                        beta_tilde,
                        intercept_tilde,
                        x_center,
                        x_scale,
                        fit_intercept,
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

          checkUserInterrupt();
      }

      if (diagnostics)
        line_searches.emplace_back(current_line_searches);

      // FISTA step
      t = 0.5*(1.0 + std::sqrt(1.0 + 4.0*t_old*t_old));
      beta = beta_tilde + (t_old - 1.0)/t * (beta_tilde - beta_tilde_old);

      if (fit_intercept)
        intercept = intercept_tilde
        + (t_old - 1.0)/t * (intercept_tilde - intercept_tilde_old);

      if (passes % 100 == 0)
        checkUserInterrupt();
    }

    double deviance = 2*family->primal(y, lin_pred);

    Results res{intercept,
                beta,
                passes,
                primals,
                duals,
                infeasibilities,
                time,
                line_searches,
                deviance};

    return res;
  }
};
