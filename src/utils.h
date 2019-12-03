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
arma::vec
linearPredictor(const arma::mat& x,
                const arma::vec& beta,
                const double intercept,
                const arma::vec& x_center,
                const arma::vec& x_scale,
                const bool standardize_features)
{
  return x*beta + intercept;
}

inline
arma::vec
linearPredictor(const arma::sp_mat& x,
                const arma::vec& beta,
                const double intercept,
                const arma::vec& x_center,
                const arma::vec& x_scale,
                const bool standardize_features)
{
  using namespace arma;

  if (standardize_features)
    return x*beta + intercept - dot(x_center/x_scale, beta);
  else
    return x*beta + intercept;
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
arma::uvec
setUnion(const arma::uvec& a,
         const arma::uvec& b)
{
  using namespace arma;

  std::vector<unsigned> out;
  std::set_union(a.begin(), a.end(),
                 b.begin(), b.end(),
                 std::back_inserter(out));

  return conv_to<uvec>::from(out);
}
