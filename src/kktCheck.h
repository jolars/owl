#pragma once

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

uvec kktCheck(const mat&   gradient,
              const mat&   beta,
              const vec&   lambda,
              const double tol)
{
  uvec nonzeros = find(beta != 0);
  uvec ord = sort_index(abs(gradient), "descend");
  vec abs_gradient_sorted = abs(gradient(ord));

  double rh = std::max(std::sqrt(datum::eps), tol*lambda(0));

  uvec out = cumsum(abs_gradient_sorted - lambda) > rh;
  out(ord) = out;
  out(nonzeros).zeros();

  umat out_mat = reshape(out, size(gradient));

  return find(any(out_mat, 1));
}
