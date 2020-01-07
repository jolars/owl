#pragma once

#include <RcppArmadillo.h>

using namespace arma;

inline void linearPredictor(mat& linear_predictor,
                            const mat& x,
                            const mat& beta,
                            const rowvec& intercept,
                            const rowvec& x_center,
                            const rowvec& x_scale,
                            const bool fit_intercept,
                            const bool standardize_features)
{
  linear_predictor = x*beta;

  if (fit_intercept) {
    for (uword k = 0; k < linear_predictor.n_cols; ++k)
      linear_predictor.col(k) += as_scalar(intercept(k));
  }
}

inline void linearPredictor(mat& linear_predictor,
                            const sp_mat& x,
                            const mat& beta,
                            const rowvec& intercept,
                            const rowvec& x_center,
                            const rowvec& x_scale,
                            const bool fit_intercept,
                            const bool standardize_features)
{
  linear_predictor = x*beta;

  for (uword k = 0; k < beta.n_cols; ++k) {
    if (fit_intercept)
      linear_predictor.col(k) += as_scalar(intercept(k));

    if (standardize_features)
      linear_predictor.col(k) -= dot(x_center/x_scale, beta.col(k));
  }
}

template <typename T>
T matrixSubset(const T& x,
               const uvec& active_set)
{
  const uword p = active_set.n_elem;
  const uword n = x.n_rows;

  T x_subset(n, p);

  for (uword j = 0; j < p; ++j) {
    uword k = active_set(j);
    x_subset.col(j) = x.col(k);
  }

  return x_subset;
}

inline uvec setUnion(const uvec& a, const uvec& b)
{
  std::vector<unsigned> out;
  std::set_union(a.begin(), a.end(),
                 b.begin(), b.end(),
                 std::back_inserter(out));

  return conv_to<uvec>::from(out);
}


inline uvec setDiff(uvec& a, uvec& b)
{
  std::vector<unsigned> out;
  std::set_difference(a.begin(), a.end(),
                      b.begin(), b.end(),
                      std::back_inserter(out));

  return conv_to<uvec>::from(out);
}

inline bool isSparse(SEXP x)
{
  bool is_sparse = false;

  if (Rf_isS4(x))
    if (Rf_inherits(x, "dgCMatrix"))
      is_sparse = true;

  return is_sparse;
}

