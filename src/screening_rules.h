#pragma once

#include <memory>
// #include <algorithm>
#include "penalties.h"
#include "families.h"

using namespace arma;

uvec
SLOPE::activeSet(const std::unique_ptr<Family>& family,
                 const mat& y,
                 const mat& gradient_prev,
                 const vec& lambda,
                 const vec& lambda_prev,
                 const std::string screening_rule)
{
  const uword p = gradient_prev.n_elem;

  uvec active_set = zeros<uvec>(p);

  if (screening_rule == "none") {

    active_set.ones();

  } else if (screening_rule == "strong") {

    vec abs_grad = abs(gradient_prev);
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
  }

  return find(active_set);
}

uvec
GroupSLOPE::activeSet(const std::unique_ptr<Family>& family,
                      const mat& y,
                      const mat& gradient_prev,
                      const vec& lambda,
                      const vec& lambda_prev,
                      const std::string screening_rule)
{
  // not implemented yet
  return regspace<uvec>(0, gradient_prev.n_elem - 1);
}
