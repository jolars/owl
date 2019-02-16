#ifndef GOLEM_FAMILIES_
#define GOLEM_FAMILIES_

#include <RcppArmadillo.h>
#include "utils.h"

class Family {
public:
  virtual
  double
  loss(const arma::vec& lin_pred) = 0;

  virtual
  arma::vec
  gradient(const arma::vec& lin_pred) = 0;

  virtual
  double
  lipschitzConstant() = 0;
};

class Gaussian : public Family {
public:
  Gaussian(const arma::mat& X, const arma::vec& y) : X(X), y(y) {}

  double
  loss(const arma::vec& lin_pred)
  {
    return 0.5*std::pow(arma::norm(lin_pred - y), 2);
  }

  arma::vec
  gradient(const arma::vec& lin_pred)
  {
    return X.t() * (lin_pred - y);
  }

  double
  lipschitzConstant()
  {
    // maximum eigenvalue of X'X
    return (arma::eig_sym(X.t() * X)).max();
  }

private:
  const arma::mat& X;
  const arma::vec& y;
};

// helper to choose family
inline
std::unique_ptr<Family>
setupFamily(const std::string& family_choice,
            const arma::mat& X,
            const arma::vec& y)
{
  return std::unique_ptr<Gaussian>(new Gaussian{X, y});
}


#endif /* GOLEM_FAMILIES_ */
