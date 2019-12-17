#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "utils.h"

using namespace arma;

class Family {
public:
  virtual
  double
  primal(const mat& y, const mat& lin_pred) = 0;

  virtual
  double
  dual(const mat& y, const mat& lin_pred) = 0;

  // this is not really the true gradient, and needs to multiplied by X'
  virtual
  mat
  pseudoGradient(const mat& y, const mat& lin_pred) = 0;

  virtual
  rowvec
  fitNullModel(const mat& y, const uword n_classes) = 0;

  std::string name = "none";
};

class Gaussian : public Family {
public:
  double
  primal(const mat& y, const mat& lin_pred)
  {
    return 0.5*pow(norm(y - lin_pred), 2);
  }

  double
  dual(const mat& y, const mat& lin_pred)
  {
    using namespace std;
    return 0.5*pow(norm(y, 2), 2) - 0.5*pow(norm(lin_pred, 2), 2);
  }

  mat
  pseudoGradient(const mat& y, const mat& lin_pred)
  {
    return -(y - lin_pred);
  }

  rowvec
  fitNullModel(const mat& y, const uword n_classes)
  {
    return mean(y);
  }

  std::string name = "gaussian";
};

class Binomial : public Family {
public:
  double
  primal(const mat& y, const mat& lin_pred)
  {
    return accu(log(1.0 + exp(-y % lin_pred)));
  }

  double
  dual(const mat& y, const mat& lin_pred)
  {
    const vec r = 1.0/(1.0 + exp(y % lin_pred));
    return as_scalar((r - 1.0).t()*log(1.0 - r) - r.t()*log(r));
  }

  mat
  pseudoGradient(const mat& y, const mat& lin_pred)
  {
    return -y / (1.0 + exp(y % lin_pred));
  }

  rowvec
  fitNullModel(const mat& y, const uword n_classes)
  {
    double pmin = 1e-9;
    double pmax = 1 - pmin;

    vec mu = clamp(mean(0.5*y + 0.5), pmin, pmax);

    return log(mu/(1 - mu));
  }

  std::string name = "binomial";
};

class Poisson : public Family {
public:
  double
  primal(const mat& y, const mat& lin_pred)
  {
    return -accu(y % lin_pred - exp(lin_pred) - lgamma(y + 1));
  }

  double
  dual(const mat& y, const mat& lin_pred)
  {
    return -accu(exp(lin_pred) % (lin_pred - 1) - lgamma(y + 1));
  }

  mat
  pseudoGradient(const mat& y, const mat& lin_pred)
  {
    return exp(lin_pred) - y;
  }

  rowvec
  fitNullModel(const mat& y, const uword n_classes)
  {
    return log(mean(y));
  }

  std::string name = "poisson";
};

class Multinomial : public Family {
public:

  double
  primal(const mat& y, const mat& lin_pred)
  {
    // logsumexp bit
    vec lp_max = max(lin_pred, 1);
    vec lse =
      log(exp(-lp_max) + sum(exp(lin_pred.each_col() - lp_max), 1)) + lp_max;

    return accu(lse) - accu(y % lin_pred);
  }

  double
  dual(const mat& y, const mat& lin_pred)
  {
    vec lp_max = max(lin_pred, 1);
    vec lse =
      log(exp(-lp_max) + sum(exp(lin_pred.each_col() - lp_max), 1)) + lp_max;

    return accu(lse) - accu(lin_pred % trunc_exp(lin_pred.each_col() - lse));
  }

  mat
  pseudoGradient(const mat& y, const mat& lin_pred)
  {
    vec lp_max = max(lin_pred, 1);
    vec lse =
      log(exp(-lp_max) + sum(exp(lin_pred.each_col() - lp_max), 1)) + lp_max;

    return exp(lin_pred.each_col() - lse) - y;
  }

  rowvec
  fitNullModel(const mat& y, const uword n_classes)
  {
    rowvec intercept = mean(y);

    return trunc_log(intercept) - accu(trunc_log(intercept))/(n_classes + 1);
  }

  std::string name = "multinomial";
};

// helper to choose family
inline
std::unique_ptr<Family>
setupFamily(const std::string& family_choice)
{
  if (family_choice == "binomial")
    return std::unique_ptr<Binomial>(new Binomial);
  else if (family_choice == "poisson")
    return std::unique_ptr<Poisson>(new Poisson);
  else if (family_choice == "multinomial")
    return std::unique_ptr<Multinomial>(new Multinomial);
  else
    return std::unique_ptr<Gaussian>(new Gaussian);
}

