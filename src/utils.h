#pragma once

#include <RcppArmadillo.h>

using namespace arma;

inline
mat
linearPredictor(const arma::mat& x,
                const arma::mat& beta,
                const rowvec& intercept,
                const arma::vec& x_center,
                const arma::vec& x_scale,
                const bool standardize_features)
{
  return (x*beta).eval().each_row() + intercept;
}

inline
mat
linearPredictor(const arma::sp_mat& x,
                const arma::mat& beta,
                const rowvec& intercept,
                const arma::vec& x_center,
                const arma::vec& x_scale,
                const bool standardize_features)
{
  uword m = beta.n_cols;

  mat lin_pred = x*beta;

  for (uword k = 0; k < m; ++k) {
    lin_pred.col(k) += intercept(k);
    if (standardize_features)
      lin_pred.col(k) -= dot(x_center/x_scale, beta.col(k));
  }

  return lin_pred;
}

template <typename T>
T
matrixSubset(const T& x,
             const arma::uvec& active_set)
{
  using namespace arma;

  const uword p = active_set.n_elem;
  const uword n = x.n_rows;

  T x_subset(n, p);

  for (uword j = 0; j < p; ++j) {
    uword k = active_set(j);
    x_subset.col(j) = x.col(k);
  }

  return x_subset;
}

inline
uvec
setUnion(const uvec& a, const uvec& b)
{
  std::vector<unsigned> out;
  std::set_union(a.begin(), a.end(),
                 b.begin(), b.end(),
                 std::back_inserter(out));

  return conv_to<uvec>::from(out);
}

inline
bool
isSparse(SEXP x)
{
  bool is_sparse = false;

  if (Rf_isS4(x))
    if (Rf_inherits(x, "dgCMatrix"))
      is_sparse = true;

  return is_sparse;
}

