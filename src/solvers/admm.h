#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "../families.h"
#include "../penalties.h"
#include "solver.h"

class ADMM : public Solver {
private:
  const bool standardize;
  const bool is_sparse;
  const bool diagnostics;
  const uword max_passes;
  const double tol_rel;
  const double tol_abs;
  const double alpha;
  const uword verbosity;

public:
  ADMM(const bool standardize,
       const bool is_sparse,
       const bool diagnostics,
       const uword max_passes,
       const double tol_rel,
       const double tol_abs,
       const double alpha,
       const uword verbosity)
       : standardize(standardize),
         is_sparse(is_sparse),
         diagnostics(diagnostics),
         max_passes(max_passes),
         tol_rel(tol_rel),
         tol_abs(tol_abs),
         alpha(alpha),
         verbosity(verbosity) {}

  Results
  fit(const sp_mat& x,
      const mat& y,
      const std::unique_ptr<Family>& family,
      const std::unique_ptr<Penalty>& penalty,
      const rowvec& intercept_init,
      const mat& beta_init,
      const bool fit_intercept,
      const vec& lambda,
      const vec& x_center,
      const vec& x_scale)
  {
    // not yet implemented
    Results res;

    Rcpp::stop("not yet implemented");

    return res;
  }

  Results
  fit(const mat& x,
      const mat& y,
      const std::unique_ptr<Family>& family,
      const std::unique_ptr<Penalty>& penalty,
      const rowvec& intercept_init,
      const mat& beta_init,
      const bool fit_intercept,
      const vec& lambda,
      const vec& x_center,
      const vec& x_scale)
  {
    using namespace arma;

    uword n = y.n_rows;
    uword p = x.n_cols;
    uword m = y.n_cols;

    rowvec intercept = intercept_init;

    mat beta(beta_init);
    mat beta_hat(beta);

    mat linear_predictor(n, m);

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

		double rho = 1/lambda.max();

		vec xTy = x.t() * y;
		mat L, U;

		if (n >= p)
		  lu(L, U, x.t()*x + rho*speye(p, p));
		else
		  lu(L, U, speye(n, n) + (1/rho)*(x * x.t()));

		vec z(p, fill::zeros);
		vec u(p, fill::zeros);
		vec z_old(z);
		double dual_feas, primal_feas;

		// ADMM loop
		uword i = 0;

		while (i < max_passes) {

		  mat q = xTy + rho*(z - u);

		  if (n >= p)
		    beta = solve(U, solve(L, q));
		  else
		    beta = q/rho - (x.t() * solve(U, solve(L, x*q)))/(rho*rho);

		  z_old = z;
		  beta_hat = alpha*beta + (1 - alpha)*z_old;
		  z = penalty->eval(beta_hat + u, lambda, 1/rho);

		  u += beta_hat - z;

			primal_feas = norm(beta - z);
			dual_feas = norm(-rho*(z - z_old));

			linear_predictor = linearPredictor(x,
                                      beta,
                                      intercept,
                                      x_center,
                                      x_scale,
                                      standardize);

			double primal =
			  family->primal(y, linear_predictor) + penalty->primal(beta, lambda);

			if (diagnostics) {
			  primals.emplace_back(primal);
        time.emplace_back(timer.toc());
        infeasibilities.emplace_back(datum::nan);
        duals.emplace_back(datum::nan);
        line_searches.emplace_back(0);
			}

			double primal_tol =
			  std::sqrt(p)*tol_abs + tol_rel*std::max(norm(beta), norm(-z));
			double dual_tol =
			  std::sqrt(p)*tol_abs + tol_rel*norm(rho*u);

			++i;

			if (primal_feas < primal_tol && dual_feas < dual_tol) {
				break;
			}
		}

		double deviance = 2*family->primal(y, linear_predictor);

    Results res{intercept,
                z,
                i,
                primals,
                duals,
                infeasibilities,
                time,
                line_searches,
                deviance};

    return res;
  }
};
