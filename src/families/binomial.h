#pragma once

#include <RcppArmadillo.h>
#include "family.h"

using namespace Rcpp;
using namespace arma;

class Binomial : public Family {
public:
  double primal(const mat& y, const mat& lin_pred)
  {
    return accu(trunc_log(1.0 + trunc_exp(-y % lin_pred)));
  }

  double dual(const mat& y, const mat& lin_pred)
  {
    const vec r = 1.0/(1.0 + trunc_exp(y % lin_pred));
    return as_scalar((r - 1.0).t()*trunc_log(1.0 - r) - r.t()*trunc_log(r));
  }

  mat pseudoGradient(const mat& y, const mat& lin_pred)
  {
    return -y / (1.0 + trunc_exp(y % lin_pred));
  }

  rowvec fitNullModel(const mat& y, const uword n_classes)
  {
    double pmin = 1e-9;
    double pmax = 1 - pmin;

    vec mu = clamp(mean(0.5*y + 0.5), pmin, pmax);

    return trunc_log(mu/(1 - mu));
  }

  std::string name()
  {
    return "binomial";
  }
};
