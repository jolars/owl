#pragma once

#include <memory>
// #include <algorithm>
#include "penalties.h"
#include "families.h"

// inline
// arma::uvec
// findInterval(const arma::vec& x,
//              const arma::vec& breaks)
// {
//   using namespace arma;
//
//   uvec out(x.n_elem);
//
//   for (uword i = 0; i < x.n_elem; ++i) {
//     auto pos = std::upper_bound(breaks.begin(), breaks.end(), x(i));
//     out[i] = std::distance(breaks.begin(), pos);
//   }
//
//   return out;
// }

arma::uvec
SLOPE::activeSet(const std::unique_ptr<Family>& family,
                 const arma::vec& y,
                 const arma::vec& gradient_prev,
                 const arma::vec& pseudo_gradient_prev,
                 const arma::vec& norms,
                 const arma::vec& lambda,
                 const arma::vec& lambda_prev,
                 const std::string screening_rule)
{
  using namespace arma;

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

  } else if (screening_rule == "safe") {

    uvec ord = sort_index(abs(gradient_prev), "descend");

    vec lh = abs(gradient_prev);

    vec rh = lambda(ord)
      - (norms*norm(pseudo_gradient_prev, 2))
        % ((lambda_prev(ord) - lambda(ord))/lambda_prev(ord));

    active_set = lh >= rh;
  }

  return find(active_set);
}

arma::uvec
GroupSLOPE::activeSet(const std::unique_ptr<Family>& family,
                      const arma::vec& y,
                      const arma::vec& gradient_prev,
                      const arma::vec& pseudo_gradient_prev,
                      const arma::vec& norms,
                      const arma::vec& lambda,
                      const arma::vec& lambda_prev,
                      const std::string screening_rule)
{
  // not implemented yet
  return arma::regspace<arma::uvec>(0, gradient_prev.n_elem - 1);
}
