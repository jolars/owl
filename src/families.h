#ifndef GOLEM_FAMILIES_
#define GOLEM_FAMILIES_

#include <RcppArmadillo.h>
#include "utils.h"

class Family {
public:
  Family(arma::mat& X, arma::vec& y) : X(X), y(y) {}

  virtual
  double
  loss(const arma::vec& lin_pred) = 0;

  // this is not really the true gradient, and needs to multiplied by X'
  virtual
  arma::vec
  gradient(const arma::vec& lin_pred) = 0;

  virtual
  double
  link(double y) = 0;

  virtual
  double
  lipschitzConstant() = 0;

  virtual
  std::pair<double, double>
  preprocessResponse(const std::string& standardize) = 0;

protected:
  arma::mat& X;
  arma::vec& y;
};

class Gaussian : public Family {
public:
  Gaussian(arma::mat& X, arma::vec& y) : Family(X, y) {}

  double
  loss(const arma::vec& lin_pred)
  {
    return 0.5*std::pow(arma::norm(lin_pred - y, 2), 2);
  }

  arma::vec
  gradient(const arma::vec& lin_pred)
  {
    return lin_pred - y;
  }

  double
  link(double y)
  {
    return y;
  }

  double
  lipschitzConstant()
  {
    // maximum eigenvalue of X'X
    return (arma::eig_sym(X.t() * X)).max();
  }

  std::pair<double, double>
  preprocessResponse(const std::string& standardize)
  {
    // always standardize gaussian responses
    double y_center = arma::mean(y);
    double y_scale = arma::stddev(y);

    y -= y_center;
    y /= y_scale;

    return std::make_pair(y_center, y_scale);
  }
};

class Binomial : public Family {
public:
  Binomial(arma::mat& X, arma::vec& y) : Family(X, y) {}

  double
  loss(const arma::vec& lin_pred)
  {
    using namespace arma;
    return accu(-y % sigmoid(lin_pred) - (1.0 - y) % sigmoid(-lin_pred));
  }

  arma::vec
  gradient(const arma::vec& lin_pred)
  {
    return sigmoid(lin_pred) - y;
  }

  double
  link(double y)
  {
    // TODO(johan): consider letting the user choose this
    double pmin = 1e-9;
    double pmax = 1.0 - pmin;
    double z = clamp(y, pmin, pmax);

    return std::log(z / (1.0 - z));
  }

  double
  lipschitzConstant()
  {
    // maximum eigenvalue of X'X
    // TODO(johan): check that this is actually true for the binomial family
    return 0.25*(arma::eig_sym(X.t() * X)).max();
  }

  std::pair<double, double>
  preprocessResponse(const std::string& standardize)
  {
    // no preprocessing for binomial response
    return std::make_pair(0.0, 1.0);
  }
};

// helper to choose family
inline
std::unique_ptr<Family>
setupFamily(const std::string& family_choice,
            arma::mat& X,
            arma::vec& y)
{
  if (family_choice == "binomial")
    return std::unique_ptr<Binomial>(new Binomial{X, y});
  else
    return std::unique_ptr<Gaussian>(new Gaussian{X, y});
}


#endif /* GOLEM_FAMILIES_ */
