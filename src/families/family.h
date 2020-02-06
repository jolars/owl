#pragma once

#include <RcppArmadillo.h>
#include "../results.h"

using namespace Rcpp;
using namespace arma;

class Family {
public:
  virtual double primal(const mat& y, const mat& lin_pred) = 0;

  virtual double dual(const mat& y, const mat& lin_pred) = 0;

  // this is not really the true gradient; it needs to multiplied by X^T
  virtual mat pseudoGradient(const mat& y, const mat& lin_pred) = 0;

  virtual rowvec fitNullModel(const mat& y, const uword n_classes) = 0;

  virtual std::string name() = 0;
};
