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
arma::vec
linearPredictor(const arma::mat& x,
                const arma::vec& beta,
                const double intercept,
                const arma::vec& x_center,
                const arma::vec& x_scale,
                const bool standardize_features,
                const bool fit_intercept)
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
                const bool standardize_features,
                const bool fit_intercept)
{
  using namespace arma;

  if (standardize_features)
    return x*beta - arma::dot(x_center/x_scale, beta) + intercept;
  else
    return x*beta + intercept;
}
