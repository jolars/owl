#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "families.h"
#include "utils.h"
#include "proxes.h"

using namespace arma;
using namespace Rcpp;

class Penalty {
public:
  virtual mat eval(const mat& beta,
                   const vec& lambda,
                   const double shrinkage) = 0;

  virtual double primal(const mat& beta, const vec& lambda) = 0;

  virtual double infeasibility(const mat& gradient, const vec& lambda) = 0;

  // returns the indices of coefficients for which the kkt test fails
  uvec kktCheck(const mat&   gradient,
                const mat&   beta,
                const vec&   lambda,
                const double tol = 1e-6)
  {
    uvec nonzeros = find(beta != 0);
    uvec ord = sort_index(abs(gradient), "descend");
    vec abs_gradient_sorted = abs(gradient(ord));

    double rh = std::max(std::sqrt(datum::eps), tol*lambda(0));

    uvec out = cumsum(abs_gradient_sorted - lambda) > rh;
    out(ord) = out;
    out(nonzeros).zeros();

    return find(out);
  }
};

class SLOPE : public Penalty {
public:
  mat eval(const mat& beta, const vec& lambda, const double shrinkage)
  {
    mat out = slopeProx(vectorise(beta), lambda, shrinkage);

    return reshape(out, size(beta));
  };

  double primal(const mat& beta, const vec& lambda)
  {
    return dot(lambda, sort(abs(vectorise(beta)), "descending"));
  }

  double infeasibility(const mat& gradient, const vec& lambda)
  {
    vec abs_gradient_sorted = sort(abs(vectorise(gradient)), "descending");
    return std::max(cumsum(abs_gradient_sorted - lambda).max(), 0.0);
  }
};

// helper to choose penalty
inline std::unique_ptr<Penalty> setupPenalty(const List& args)
{
  return std::unique_ptr<SLOPE>(new SLOPE);
}
