#pragma once

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

template <typename T>
vec lambdaMax(const T& x,
              const mat& y,
              const rowvec& x_center,
              const rowvec& x_scale,
              const rowvec& y_scale,
              const uword n_targets,
              const std::string& family,
              const bool standardize_features,
              const bool is_sparse)
{
  const uword p = x.n_cols;
  mat lambda_max(p, n_targets);

  if (family == "binomial") {
    vec y_new = (y + 1)/2;

    // standardize
    double y_center = mean(y_new);
    y_new -= y_center;

    lambda_max = x.t() * y_new;

    if (is_sparse && standardize_features) {
      for (uword j = 0; j < p; ++j)
        lambda_max(j) -= accu(y_new * x_center(j)/x_scale(j));
    }

  } else if (family == "multinomial") {

    rowvec y_bar = mean(y);
    rowvec y_std = stddev(y, 1);

    mat y_map = y;

    for (uword k = 0; k < n_targets; ++k) {
      y_map.col(k) -= y_bar(k);
      y_map.col(k) /= y_std(k);
    }

    lambda_max = x.t() * y_map;

    for (uword k = 0; k < n_targets; ++k) {
      lambda_max.col(k) *= y_std(k);
    }

  } else if (family == "poisson") {

    lambda_max = x.t() * (1 - y);

    if (is_sparse && standardize_features) {
      for (uword j = 0; j < p; ++j)
        lambda_max(j) -= accu((1 - y) * x_center(j)/x_scale(j));
    }

  } else {

    lambda_max = x.t() * y;

    if (is_sparse && standardize_features) {
      for (uword j = 0; j < p; ++j)
        lambda_max(j) -= accu(y * x_center(j)/x_scale(j));
    }
  }

  return abs(vectorise(lambda_max));
}
