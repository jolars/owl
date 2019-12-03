#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "utils.h"

class Family {
protected:
  const bool fit_intercept;
  const bool standardize;

public:
  Family(const bool fit_intercept, const bool standardize)
         : fit_intercept(fit_intercept), standardize(standardize) {}

  virtual
  double
  primal(const arma::vec& y, const arma::vec& lin_pred) = 0;

  virtual
  double
  dual(const arma::vec& y, const arma::vec& lin_pred) = 0;

  // this is not really the true gradient, and needs to multiplied by X'
  virtual
  arma::vec
  pseudoGradient(const arma::vec& y, const arma::vec& lin_pred) = 0;

  virtual
  double
  link(const double y) = 0;

  virtual
  double
  fitNullModel(const arma::vec& y) = 0;
};

class Gaussian : public Family {
public:
  Gaussian(const bool fit_intercept, const bool standardize)
           : Family(fit_intercept, standardize) {}

  double
  primal(const arma::vec& y, const arma::vec& lin_pred)
  {
    return 0.5*std::pow(arma::norm(y - lin_pred), 2);
  }

  double
  dual(const arma::vec& y, const arma::vec& lin_pred)
  {
    using namespace arma;
    using namespace std;
    return 0.5*pow(norm(y, 2), 2) - 0.5*pow(norm(lin_pred, 2), 2);
  }

  arma::vec
  pseudoGradient(const arma::vec& y, const arma::vec& lin_pred)
  {
    return -(y - lin_pred);
  }

  double
  link(const double y)
  {
    return y;
  }

  double
  fitNullModel(const arma::vec& y)
  {
    return arma::mean(y);
  }
};

class Binomial : public Family {
public:
  Binomial(const bool fit_intercept, const bool standardize)
           : Family(fit_intercept, standardize) {}

  double
  primal(const arma::vec& y, const arma::vec& lin_pred)
  {
    using namespace arma;
    return accu(log(1.0 + exp(-y % lin_pred)));
  }

  double
  dual(const arma::vec& y, const arma::vec& lin_pred)
  {
    using namespace arma;
    const arma::vec r = 1.0/(1.0 + arma::exp(y % lin_pred));
    return arma::as_scalar((r - 1.0).t()*log(1.0 - r) - r.t()*log(r));
  }

  arma::vec
  pseudoGradient(const arma::vec& y, const arma::vec& lin_pred)
  {
    return -y / (1.0 + arma::exp(y % lin_pred));
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

  double
  fitNullModel(const arma::vec& y)
  {
    using namespace arma;

    return std::log(mean(y*0.5 + 0.5)/(1 - mean(y*0.5 + 0.5)));
  }
};

// helper to choose family
inline
std::unique_ptr<Family>
setupFamily(const std::string& family_choice,
            const bool fit_intercept,
            const bool standardize)
{
  if (family_choice == "binomial")
    return std::unique_ptr<Binomial>(new Binomial{fit_intercept,
                                                  standardize});
  else
    return std::unique_ptr<Gaussian>(new Gaussian{fit_intercept,
                                                  standardize});
}

