#pragma once

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

inline double Infeasibility(const mat& gradient, const vec& lambda)
{
  return std::max(cumsum(sort(abs(vectorise(gradient)),
                              "descending") - lambda).max(), 0.0);
}
