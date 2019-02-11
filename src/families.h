#ifndef GOLEM_FAMILIES_
#define GOLEM_FAMILIES_

#include <RcppArmadillo.h>

class Family {
public:
  virtual
  double
  loss(const arma::vec& lin_pred) = 0;

  virtual
  arma::vec
  gradient(const arma::vec& lin_pred) = 0;
};

class Gaussian : public Family {
public:
  Gaussian(const arma::mat& X, const arma::vec& y) : X(X), y(y) {}

  double
  loss(const arma::vec& lin_pred)
  {
    const arma::vec res = lin_pred - y;
    return 0.5*arma::dot(res, res);
  }

  arma::vec
  gradient(const arma::vec& lin_pred)
  {
    return X.t() * (lin_pred - y);
  }

private:
  const arma::mat& X;
  const arma::vec& y;
};


#endif /* GOLEM_FAMILIES_ */
