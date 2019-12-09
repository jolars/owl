#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "utils.h"

class Family {
public:
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
  fitNullModel(const arma::vec& y) = 0;
};

class Gaussian : public Family {
public:
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
  fitNullModel(const arma::vec& y)
  {
    return arma::mean(y);
  }
};

class Binomial : public Family {
public:
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
  fitNullModel(const arma::vec& y)
  {
    using namespace arma;

    return std::log(mean(y*0.5 + 0.5)/(1 - mean(y*0.5 + 0.5)));
  }
};

class Poisson : public Family {
public:
  double
  primal(const arma::vec& y, const arma::vec& lin_pred)
  {
    using namespace arma;
    return -accu(y % lin_pred - exp(lin_pred));
  }

  double
  dual(const arma::vec& y, const arma::vec& lin_pred)
  {
    using namespace arma;
    const vec theta = y - exp(lin_pred);
    return -accu((y - theta) % log(y - theta) - y + theta);
  }

  arma::vec
  pseudoGradient(const arma::vec& y, const arma::vec& lin_pred)
  {
    return arma::exp(lin_pred) - y;
  }

  double
  fitNullModel(const arma::vec& y)
  {
    return std::log(arma::mean(y));
  }
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
  else
    return std::unique_ptr<Gaussian>(new Gaussian);
}

