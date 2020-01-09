#pragma once

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

void standardize(mat& x, rowvec& x_center, rowvec& x_scale)
{
  const uword p = x.n_cols;

  for (uword j = 0; j < p; ++j) {
    x_center(j) = mean(x.col(j));
    x_scale(j) = norm(x.col(j) - x_center(j));
  }

  // don't scale zero-variance predictors
  x_scale.replace(0, 1);

  for (uword j = 0; j < p; ++j) {
    x.col(j) -= x_center(j);
    x.col(j) /= x_scale(j);
  }
}

void standardize(sp_mat& x, rowvec& x_center, rowvec& x_scale)
{
  const uword p = x.n_cols;
  const uword n = x.n_rows;

  for (uword j = 0; j < p; ++j) {
    x_center(j) = accu(x.col(j))/n;
    x_scale(j) = norm(x.col(j) - x_center(j));
  }

  // don't scale zero-variance predictors
  x_scale.replace(0, 1);

  for (uword j = 0; j < p; ++j)
    x.col(j) /= x_scale(j);
}
