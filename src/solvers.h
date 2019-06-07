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
          std::vector<double> primals,
          std::vector<double> duals,
          std::vector<double> infeasibilities,
          std::vector<double> time)
          : intercept(intercept),
            beta(beta),
            passes(passes),
            primals(primals),
            duals(duals),
            infeasibilities(infeasibilities),
            time(time) {}

  arma::rowvec intercept;
  arma::mat beta;
  arma::uword passes;
  std::vector<double> primals;
  std::vector<double> duals;
  std::vector<double> infeasibilities;
  std::vector<double> time;
};

class Solver {
protected:
  bool diagnostics;
  std::vector<double> primals;
  std::vector<double> duals;
  std::vector<double> infeasibilities;
  std::vector<double> time;
  arma::uword path_iter = 0;
  arma::rowvec intercept;
  arma::mat beta;
};


class FISTA : public Solver {
private:
  double L;
  const bool standardize;
  arma::vec x_scaled_center;
  const bool is_sparse;
  arma::uword max_passes;
  const double eta = 2.0;
  double tol_rel_gap = 1e-6;
  double tol_infeas = 1e-6;

public:
  FISTA(arma::rowvec&& intercept_init,
        arma::mat&& beta_init,
        const double lipschitz_constant,
        const bool standardize,
        const arma::vec x_scaled_center,
        const bool is_sparse,
        const Rcpp::List& args)
        : L(lipschitz_constant),
          standardize(standardize),
          x_scaled_center(x_scaled_center),
          is_sparse(is_sparse)
  {
    using Rcpp::as;

    intercept   = std::move(intercept_init);
    beta        = std::move(beta_init);

    max_passes  = as<arma::uword>(args["max_passes"]);
    diagnostics = as<bool>(args["diagnostics"]);
    tol_rel_gap = as<double>(args["tol_rel_gap"]);
    tol_infeas  = as<double>(args["tol_infeas"]);
  }

  template <typename T>
  Results
  fit(const T& x,
      const arma::mat& y,
      const std::unique_ptr<Family>& family,
      const std::unique_ptr<Penalty>& penalty,
      const bool fit_intercept)
  {
    using namespace arma;

    uword n = x.n_rows;
    uword p = x.n_cols;
    uword m = y.n_cols;

    rowvec intercept_tilde(intercept);
    rowvec intercept_tilde_old(intercept_tilde);

    mat beta_tilde(beta);
    mat beta_tilde_old(beta_tilde);

    mat lin_pred(n, m);

    mat g(p, m, fill::zeros);
    mat pseudo_g(g);
    rowvec g_intercept(m, fill::zeros);

    double t = 1;

    uword i = 0;
    bool accepted = false;

    tol_infeas *= penalty->lambdaInfeas();
    tol_infeas =
      tol_infeas < std::sqrt(datum::eps) ? std::sqrt(datum::eps) : tol_infeas;

    // diagnostics
    wall_clock timer;

    if (diagnostics) {
      primals.reserve(max_passes);
      duals.reserve(max_passes);
      infeasibilities.reserve(max_passes);
      time.reserve(max_passes);
      timer.tic();
    }

    family->eval(x, y, intercept, beta, x_scaled_center);

    while (!accepted && i < max_passes) {
      // gradient
      double f = family->primal();
      pseudo_g = family->gradient(y);
      g = x.t() * pseudo_g;

      // adjust gradient if sparse and standardizing
      if (is_sparse && standardize) {
        for (uword j = 0; j < p; ++j)
          g(j) -= accu(x_scaled_center(j)*pseudo_g);
      }

      if (fit_intercept) {
        g_intercept = mean(pseudo_g);
        intercept_tilde_old = intercept_tilde;
      }

      double primal = f + penalty->primal(beta_tilde);
      double dual = family->dual(y);
      double infeasibility = penalty->infeasibility(g);

      accepted = (std::abs(primal - dual)/std::max(1.0, primal) < tol_rel_gap)
                  && (infeasibility <= tol_infeas);

      if (diagnostics) {
        time.push_back(timer.toc());
        primals.push_back(primal);
        infeasibilities.push_back(infeasibility);
        duals.push_back(dual);
      }

      if (accepted)
        break;

      beta_tilde_old = beta_tilde;
      double f_old = f;
      double t_old = t;
      beta_tilde_old = beta_tilde;

      // // Lipschitz search
      while (true) {
        // Update beta and intercept
        beta_tilde = penalty->eval(beta - (1.0/L)*g, 1.0/L);

        mat d = beta_tilde - beta;
        if (fit_intercept)
          intercept_tilde = intercept - (1.0/L)*g_intercept;

        family->eval(x, y, intercept_tilde, beta_tilde, x_scaled_center);

        f = family->primal();
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

      family->eval(x, y, intercept, beta, x_scaled_center);

      i++;
    }

    Results res{intercept, beta, i, primals, duals, infeasibilities, time};

    return res;
  }
};

