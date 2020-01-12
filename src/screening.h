#pragma once

#include <RcppArmadillo.h>

using namespace arma;

uvec activeSet(const mat& gradient_prev,
               const vec& lambda,
               const vec& lambda_prev)
{
  const uword p = gradient_prev.n_elem;

  uvec active_set = zeros<uvec>(p);

  vec abs_grad = abs(vectorise(gradient_prev));
  uvec ord = sort_index(abs_grad, "descend");
  vec abs_grad_sorted = abs_grad(ord);
  vec tmp = abs_grad_sorted + lambda_prev - 2*lambda;

  uword i = 0;
  uword k = 0;

  double s = 0;

  while (i + k < p) {
    s += tmp(k+i);

    if (s >= 0) {
      k = k + i + 1;
      i = 0;
      s = 0;
    } else {
      i++;
    }
  }

  active_set.head(k).ones();
  active_set(ord) = active_set;

  umat active_set_mat = reshape(active_set, size(gradient_prev));

  return find(any(active_set_mat, 1));
}
