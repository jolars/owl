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
  arma::uword max_passes = 1e3;
  double tol_rel_gap = 1e-6;
  double tol_infeas = 1e-6;

public:
  ADMM(const bool standardize,
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

  Results
  fit(const arma::sp_mat& x,
      const arma::vec& y,
      const std::unique_ptr<Family>& family,
      const std::unique_ptr<Penalty>& penalty,
      const double intercept_init,
      const arma::vec& beta_init,
      const bool fit_intercept,
      const double lipschitz_constant,
      const arma::vec& lambda,
      const arma::vec& x_center,
      const arma::vec& x_scale)
  {
    // not yet implemented
    Results res;

    return res;
  }

  Results
  fit(const arma::mat& x,
      const arma::vec& y,
      const std::unique_ptr<Family>& family,
      const std::unique_ptr<Penalty>& penalty,
      const double intercept_init,
      const arma::vec& beta_init,
      const bool fit_intercept,
      const double lipschitz_constant,
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
    // vec beta_tilde(beta);
    // vec beta_tilde_old(beta_tilde);

    // vec lin_pred(n);

    // vec gradient(p, fill::zeros);
    // vec pseudo_gradient(n, fill::zeros);
    // double gradient_intercept = 0.0;

    double learning_rate = 1.0/lipschitz_constant;
    // double learning_rate = 1;

    // line search parameters
    // double eta = 0.5;

    // FISTA parameters
    // double t = 1;

    // uword passes = 0;
    // bool accepted = false;

    // diagnostics
    // wall_clock timer;
    // std::vector<double> primals;
    // std::vector<double> duals;
    // std::vector<double> infeasibilities;
    // std::vector<double> time;
    // std::vector<unsigned> line_searches;
    //
    // if (diagnostics) {
    //   primals.reserve(max_passes);
    //   duals.reserve(max_passes);
    //   infeasibilities.reserve(max_passes);
    //   time.reserve(max_passes);
    //   line_searches.reserve(max_passes);
    //   timer.tic();
    // }

		double rho = 1/lambda.max();
		double tol_inf = 1e-8;

    // Precompute M = (X^TX + rho I)^{-1}
		// and MXtY = M * X^T * Y for proximal steps of quadratic part
		mat M;
		M = x.t() * x;
		M.diag() += rho;
		M = inv(M);
		vec MXtY = M * (x.t() * y);
		vec lam_seq_rho = lambda/rho;

		// Prepare variables before starting ADMM loop
		vec z = zeros(p);
		vec u = zeros(p);
		vec z_new(p);
		double dual_feas, primal_feas;

		// ADMM loop
		uword i = 0;

		for (i = 0; i < max_passes; ++i) {

			beta = MXtY + M*rho*(z - u);
			z_new = penalty->eval(beta + u, lambda, 1/rho);
			u += beta - z_new;

			dual_feas = norm(rho*(z_new - z), 2);
			primal_feas = norm(z_new - beta, 2);

			z = z_new;

			if (primal_feas < tol_inf && dual_feas < tol_inf) {
				break;
			}
		}

    Results res{intercept,
                beta,
                i + 1,
                primals,
                duals,
                infeasibilities,
                time,
                line_searches};

    return res;
  }
};
//
//
// using namespace arma;
// using namespace Rcpp;
//
// arma::vec slope_admm(const mat& X,
//                      const vec& Y,
//                      NumericVector& lambda,
//                      const int& p,
//                      const double& rho,
//                      int max_iter=500,
//                      double tol_inf = 1e-08) {
//
// 		// Precompute M = (X^TX + rho I)^{-1}
// 		// and MXtY = M * X^T * Y for proximal steps of quadratic part
// 		mat M = X.t() * X;
// 		for (int i=0; i<p; ++i) {
// 			M.at(i,i) += rho;
// 		}
// 		M = M.i();
// 		vec MXtY = M * (X.t() * Y);
// 		NumericVector lam_seq_rho = lambda/rho;
//
// 		// Prepare variables before starting ADMM loop
// 		int i=0;
// 		vec x = zeros(p);
// 		vec z = zeros(p);
// 		vec u = zeros(p);
// 		NumericVector z_new = NumericVector(p);
// 		vec z_new_arma = zeros(p);
// 		NumericVector x_plus_u(p);
// 		double dual_feas, primal_feas;
//
// 		// ADMM loop
// 		while (i < max_iter) {
// 			x = MXtY + M*(rho*(z - u));
// 			x_plus_u = as<NumericVector>(wrap(x+u));
// 			z_new = prox_sorted_L1_C(x_plus_u, lam_seq_rho);
// 			z_new_arma = as<arma::vec>(z_new);
// 			u += (x - z_new_arma);
//
// 			dual_feas = arma::norm(rho*(z_new_arma - z));
// 			primal_feas = arma::norm(z_new_arma - x);
//
// 			z = z_new_arma;
// 			if (primal_feas < tol_inf && dual_feas < tol_inf){
// 				i = max_iter;
// 			}
//
// 			++i;
// 		}
//
// 		return z;
// }
