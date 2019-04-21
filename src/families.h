#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "utils.h"

class Family {
public:
  virtual
  double
  primal(const arma::mat& lin_pred, const arma::mat& y) = 0;

  virtual
  double
  dual(const arma::mat& lin_pred, const arma::mat&y) = 0;

  // this is not really the true gradient, and needs to multiplied by X'
  virtual
  arma::mat
  gradient(const arma::mat& lin_pred, const arma::mat& y) = 0;

  virtual
  double
  link(const double y) = 0;
};

class Gaussian : public Family {
public:
  Gaussian() {};

  double
  primal(const arma::mat& lin_pred, const arma::mat& y)
  {
    return 0.5*std::pow(arma::norm(lin_pred - y), 2);
  }

  double
  dual(const arma::mat& lin_pred, const arma::mat& y)
  {
    return -primal(lin_pred, y) - arma::dot(lin_pred - y, y);
  }

  arma::mat
  gradient(const arma::mat& lin_pred, const arma::mat& y)
  {
    return lin_pred - y;
  }

  double
  link(const double y)
  {
    return y;
  }
};

class Binomial : public Family {
public:
  Binomial() {};

  double
  primal(const arma::mat& lin_pred, const arma::mat& y)
  {
    return arma::accu(arma::log(1.0 + arma::exp(-y % lin_pred)));
  }

  double
  dual(const arma::mat& lin_pred, const arma::mat&y)
  {
    using namespace arma;
    const arma::vec r = 1.0/(1.0 + arma::exp(y % lin_pred));
    return arma::as_scalar((r-1.0).t()*log(1.0-r) - r.t()*log(r));
  }

  arma::mat
  gradient(const arma::mat& lin_pred, const arma::mat& y)
  {
    return -y / (arma::exp(y % lin_pred) + 1.0);
  }

  double
  link(const double y)
  {
    // TODO(johan): consider letting the user choose this
    double pmin = 1e-9;
    double pmax = 1.0 - pmin;
    double z = clamp(y, pmin, pmax);

    return std::log((y + 1.0)/2 / (1.0 - z));
  }
};

// helper to choose family
inline
std::unique_ptr<Family>
setupFamily(const std::string& family_choice)
{
  if (family_choice == "binomial")
    return std::unique_ptr<Binomial>(new Binomial{});
  else
    return std::unique_ptr<Gaussian>(new Gaussian{});
}

