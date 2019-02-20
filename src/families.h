#ifndef GOLEM_FAMILIES_
#define GOLEM_FAMILIES_

#include <RcppArmadillo.h>
#include "utils.h"

class Family {
public:
  virtual
  double
  loss(const arma::vec& lin_pred, const arma::vec& y) = 0;

  // this is not really the true gradient, and needs to multiplied by X'
  virtual
  arma::vec
  gradient(const arma::vec& lin_pred, const arma::vec& y) = 0;

  virtual
  double
  link(const double y) = 0;

  virtual
  double
  lipschitzConstant(const arma::mat& X) = 0;

  virtual
  void
  preprocessResponse(arma::vec& y,
                     double& y_center,
                     double& y_scale,
                     const std::string& standardize) = 0;
};

class Gaussian : public Family {
public:
  Gaussian() {};

  double
  loss(const arma::vec& lin_pred, const arma::vec& y)
  {
    const arma::vec residual = lin_pred - y;
    return 0.5*arma::dot(residual, residual);
  }

  arma::vec
  gradient(const arma::vec& lin_pred, const arma::vec& y)
  {
    return lin_pred - y;
  }

  double
  link(const double y)
  {
    return y;
  }

  double
  lipschitzConstant(const arma::mat& X)
  {
    // maximum eigenvalue of X'X
    // TODO(johan): account for the intercept
    return (arma::eig_sym(X.t() * X)).max();
  }

  void
  preprocessResponse(arma::vec& y,
                     double& y_center,
                     double& y_scale,
                     const std::string& standardize)
  {
    // always standardize gaussian responses
    y_center = arma::mean(y);
    y_scale = arma::stddev(y);

    y -= y_center;
    y /= y_scale;
  }
};

class Binomial : public Family {
public:
  Binomial() {};

  double
  loss(const arma::vec& lin_pred, const arma::vec& y)
  {
    using namespace arma;
    return accu(-y % sigmoid(lin_pred) - (1.0 - y) % sigmoid(-lin_pred));
  }

  arma::vec
  gradient(const arma::vec& lin_pred, const arma::vec& y)
  {
    return sigmoid(lin_pred) - y;
  }

  double
  link(const double y)
  {
    // TODO(johan): consider letting the user choose this
    double pmin = 1e-9;
    double pmax = 1.0 - pmin;
    double z = clamp(y, pmin, pmax);

    return std::log(z / (1.0 - z));
  }

  double
  lipschitzConstant(const arma::mat& X)
  {
    // maximum eigenvalue of X'X
    // TODO(johan): check that this is actually true for the binomial family
    return 0.25*(arma::eig_sym(X.t() * X)).max();
  }

  void
  preprocessResponse(arma::vec& y,
                     double& y_center,
                     double& y_scale,
                     const std::string& standardize)
  {
    // no preprocessing for binomial response
    y_center = 0;
    y_scale = 1;
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


#endif /* GOLEM_FAMILIES_ */
