#pragma once

#include <RcppArmadillo.h>

using namespace arma;

uvec activeSet(const mat& y,
               const mat& gradient_prev,
               const vec& lambda,
               const vec& lambda_prev)
{
  const uword p = gradient_prev.n_elem;

  uvec active_set = zeros<uvec>(p);

  vec abs_grad = abs(vectorise(gradient_prev));
  uvec ord = sort_index(abs_grad, "descend");
  vec abs_grad_sorted = abs_grad(ord);

  uword i = 0;
  uword k = 0;

  double s{0};

  while (i + k < p) {

    s += abs_grad_sorted(k) + lambda_prev(k) - 2*lambda(k);

    if (s >= 0) {
      k = k + i + 1;
      i = 0;
      s = 0;
    } else {
      i++;
    }
  }

  active_set.head(k).fill(1);
  active_set(ord) = active_set;

  umat active_set_mat = reshape(active_set, size(gradient_prev));

  return find(any(active_set_mat, 1));
}
