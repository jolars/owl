#pragma once

#include <RcppArmadillo.h>

inline
arma::vec
sigmoid(const arma::vec& x)
{
  return 1.0/(1.0 + arma::exp(-x));
}

inline
double
sigmoid(const double x)
{
  return 1.0/(1.0 + std::exp(-x));
}

//' Clamp a value to [min, max]
//' @param x value to clamp
//' @param min min
//' @param max max
//' @noRd
template <typename T>
inline
T
clamp(const T& x, const T& min, const T& max)
{
  return x > max ? max : (x < min ? min : x);
}

inline
arma::mat
innerProduct(const arma::mat& x,
             const arma::mat& beta,
             const arma::vec& x_scaled_center)
{
  return x * beta;
}

inline
arma::mat
innerProduct(const arma::sp_mat& x,
             const arma::mat& beta,
             const arma::vec& x_scaled_center)
{
  using namespace arma;

  uword n = x.n_rows;
  uword p = x.n_cols;
  uword m = beta.n_cols;

  mat lin_pred(m, n, fill::zeros);

  sp_mat::const_col_iterator it, it_end;

  for (uword j = 0; j < p; ++j) {
    it = x.begin_col(j);
    it_end = x.end_col(j);

    for (; it != it_end; ++it) {
      uword i = it.row();

      for (uword k = 0; k < m; ++k) {
        lin_pred(k, i) += ((*it) - x_scaled_center(j))*beta(j, k);
      }
    }
  }

  return lin_pred.t();
}

// class ConvergenceCheck {
// public:
//   ConvergenceCheck(const arma::rowvec& intercept_old,
//                    const arma::vec& beta_old,
//                    const double tol)
//                    : intercept_old(intercept_old),
//                      beta_old(beta_old),
//                      tol(tol) {}
//
//   bool
//   operator()(const arma::rowvec& intercept_new,
//              const arma::vec& beta_new)
//   {
//     double max_change = std::max(abs(intercept_new - intercept_old).max(),
//                                  abs(beta_new - beta_old).max());
//     double max_size   = std::max(abs(intercept_new).max(),
//                                  abs(beta_new).max());
//
//     bool all_zero  = (max_size == 0.0) && (max_change == 0.0);
//     bool no_change = (max_size != 0.0) && (max_change/max_size <= tol);
//
//     beta_old = beta_new;
//
//     return all_zero || no_change;
//   }
//
// private:
//   arma::rowvec intercept_old;
//   arma::vec beta_old;
//   const double tol;
// };
