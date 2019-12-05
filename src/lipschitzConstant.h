#pragma once

#include <RcppArmadillo.h>

double
lipschitzConstant(const arma::mat& x,
                  const arma::vec& x_center,
                  const arma::vec& x_scale,
                  const bool standardize_features,
                  const std::string family)
{
  if (x.n_rows > x.n_cols) {
    return arma::eig_sym(x.t() * x).max();
  } else {
    return arma::eig_sym(x * x.t()).max();
  }
}

double
lipschitzConstant(const arma::sp_mat& x,
                  const arma::vec& x_center,
                  const arma::vec& x_scale,
                  const bool standardize_features,
                  const std::string family)
{
  using namespace arma;

  const uword p = x.n_cols;
  const uword n = x.n_rows;

  double out{0};

  if (standardize_features) {
    // x is scaled but not centered, so we need to compute
    // x^T x - x^T x_cs - x_cs^T x + x_cs^T x_cs (for n > p) where x_cs is
    // a matrix where each column is the scaled centers (means) of columns

    mat xx;

    const vec x_center_scaled = x_center/x_scale;

    if (n > p) {

      xx = x.t() * x;

      for (uword i = 0; i < p; ++i) {
        for (uword j = 0; j < p; ++j) {
          xx(i, j) -= 2*accu(x_center_scaled(i)*x.col(j));
          xx(i, j) += x_center_scaled(i)*x_center_scaled(j)*n;
        }
      }
    } else {

      xx = x * x.t();

      for (uword i = 0; i < n; ++i) {
        for (uword j = 0; j < n; ++j) {
          // NOTE(JL): this is not efficient right now. Can we get rid of
          // x.row(j) somehow?
          double tmp = dot(x_center_scaled.t(), x.row(j));
          xx(i, j) -= tmp;
          xx(j, i) -= tmp;
        }
      }
      xx += accu(square(x_center_scaled));
    }

    out = eig_sym(xx).max();

  } else {

    sp_mat xx;

    if (n > p)
      xx = x.t() * x;
    else
      xx = x * x.t();

    out = as_scalar(eigs_sym(xx, 1));
  }

  if (family == "binomial")
    out /= 4;

  return out;
}
