#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "utils.h"

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

class Gaussian : public Family {
public:
  double primal(const mat& y, const mat& lin_pred)
  {
    return 0.5*pow(norm(y - lin_pred), 2);
  }

  double dual(const mat& y, const mat& lin_pred)
  {
    using namespace std;
    return 0.5*pow(norm(y, 2), 2) - 0.5*pow(norm(lin_pred, 2), 2);
  }

  mat pseudoGradient(const mat& y, const mat& lin_pred)
  {
    return lin_pred - y;
  }

  rowvec fitNullModel(const mat& y, const uword n_classes)
  {
    return mean(y);
  }

  std::string name()
  {
    return "gaussian";
  }
};

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

class Poisson : public Family {
public:
  double primal(const mat& y, const mat& lin_pred)
  {
    return -accu(y % lin_pred - trunc_exp(lin_pred) - lgamma(y + 1));
  }

  double dual(const mat& y, const mat& lin_pred)
  {
    return -accu(trunc_exp(lin_pred) % (lin_pred - 1) - lgamma(y + 1));
  }

  mat pseudoGradient(const mat& y, const mat& lin_pred)
  {
    return trunc_exp(lin_pred) - y;
  }

  rowvec fitNullModel(const mat& y, const uword n_classes)
  {
    return trunc_log(mean(y));
  }

  std::string name()
  {
    return "poisson";
  }
};

class Multinomial : public Family {
public:

  double primal(const mat& y, const mat& lin_pred)
  {
    // logsumexp bit
    vec lp_max = max(lin_pred, 1);
    vec lse =
      trunc_log(exp(-lp_max) + sum(trunc_exp(lin_pred.each_col() - lp_max), 1)) + lp_max;

    return accu(lse) - accu(y % lin_pred);
  }

  double dual(const mat& y, const mat& lin_pred)
  {
    vec lp_max = max(lin_pred, 1);
    vec lse =
      trunc_log(exp(-lp_max) + sum(trunc_exp(lin_pred.each_col() - lp_max), 1)) + lp_max;

    return accu(lse) - accu(lin_pred % trunc_exp(lin_pred.each_col() - lse));
  }

  mat pseudoGradient(const mat& y, const mat& lin_pred)
  {
    vec lp_max = max(lin_pred, 1);
    vec lse =
      trunc_log(exp(-lp_max) + sum(trunc_exp(lin_pred.each_col() - lp_max), 1)) + lp_max;

    return trunc_exp(lin_pred.each_col() - lse) - y;
  }

  rowvec fitNullModel(const mat& y, const uword n_classes)
  {
    const uword m = y.n_cols;

    rowvec mu = mean(y);
    rowvec log_mu = trunc_log(mu);

    return log_mu - accu(log_mu + trunc_log(1 - accu(mu)))/(m + 1);
  }

  std::string name()
  {
    return "multinomial";
  }
};

// helper to choose family
inline std::unique_ptr<Family> setupFamily(const std::string& family_choice)
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

