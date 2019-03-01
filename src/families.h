#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "utils.h"

class Family {
public:
  virtual
  double
  loss(const arma::mat& lin_pred, const arma::mat& y) = 0;

  // this is not really the true gradient, and needs to multiplied by X'
  virtual
  arma::mat
  gradient(const arma::mat& lin_pred, const arma::mat& y) = 0;

  virtual
  double
  link(const double y) = 0;

  virtual
  double
  lipschitzConstant(const arma::mat& x, const bool fit_intercept) = 0;

  virtual
  void
  preprocessResponse(arma::mat& y,
                     arma::rowvec& y_center,
                     arma::rowvec& y_scale,
                     const std::string& standardize) = 0;

  virtual
  double
  lambdaMax(const arma::mat& x,
            const arma::mat& y,
            const arma::rowvec& y_scale) = 0;

};

class Gaussian : public Family {
public:
  Gaussian() {};

  double
  loss(const arma::mat& lin_pred, const arma::mat& y)
  {
    return 0.5*std::pow(arma::norm(lin_pred - y), 2);
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

  double
  lipschitzConstant(const arma::mat& x, const bool fit_intercept)
  {
    using namespace arma;

    double row_squared_l2_norm_max = sum(square(x), 1).max();

    return row_squared_l2_norm_max + static_cast<double>(fit_intercept);
  }

  void
  preprocessResponse(arma::mat& y,
                     arma::rowvec& y_center,
                     arma::rowvec& y_scale,
                     const std::string& standardize)
  {
    // always center gaussian responses
    y_center = arma::mean(y);
    // y_scale = arma::stddev(y);
    y_scale.ones();

    y -= y_center(0);
    y /= y_scale(0);
  }

  double
  lambdaMax(const arma::mat& x,
            const arma::mat& y,
            const arma::rowvec& y_scale)
  {
    using namespace arma;
    return y_scale(0)*abs(x.t()*y.col(0)).max()/y.n_rows;
  }
};

class Binomial : public Family {
public:
  Binomial() {};

  double
  loss(const arma::mat& lin_pred, const arma::mat& y)
  {
    return arma::accu(arma::log(1.0 + arma::exp(-y % lin_pred)));
  }

  arma::mat
  gradient(const arma::mat& lin_pred, const arma::mat& y)
  {
    return -y/(arma::exp(y%lin_pred) + 1.0);
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
  lipschitzConstant(const arma::mat& x, const bool fit_intercept)
  {
    using namespace arma;

    double row_squared_l2_norm_max = sum(square(x), 1).max();

    return 0.25*(row_squared_l2_norm_max + static_cast<double>(fit_intercept));
  }

  void
  preprocessResponse(arma::mat& y,
                     arma::rowvec& y_center,
                     arma::rowvec& y_scale,
                     const std::string& standardize)
  {
    // no preprocessing for binomial response
    y_center.zeros();
    y_scale.ones();
  }

  double
  lambdaMax(const arma::mat& x,
            const arma::mat& y,
            const arma::rowvec& y_scale)
  {
    using namespace arma;

    auto n = y.n_rows;

    // convert y that is is {-1, -1} to {0, 1}
    vec y_tmp = (vectorise(y)+ 1.0)/2.0;

    double y_bar = mean(y_tmp);
    double y_sd = stddev(y_tmp);

    y_tmp -= y_bar;
    y_tmp /= y_sd != 0.0 ? y_sd : 1.0;

    return y_sd*abs(x.t() * y_tmp).max()/n;
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

