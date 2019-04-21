#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "families.h"
#include "utils.h"

class Penalty {
public:
  virtual
  arma::mat
  eval(const arma::mat& y, const double L) = 0;

  virtual
  double
  primal(const arma::mat& beta) = 0;

  virtual
  double
  infeasibility(const arma::mat& grad) = 0;

  virtual
  void
  step(const arma::uword i) = 0;

  virtual
  arma::uword
  pathLength() = 0;

  virtual
  double
  lambdaInfeas() = 0;
};

class SLOPE : public Penalty {
public:
  const arma::vec lambda;
  const double sigma;

  SLOPE(const double sigma, const arma::vec lambda)
        : sigma(sigma), lambda(lambda) {}

  arma::mat
  eval(const arma::mat& beta, const double step_size)
  {
    using namespace arma;

    uword p = beta.n_rows;

    // collect sign of beta and work with sorted absolutes
    mat beta_sign = sign(beta);
    vec beta2 = abs(beta);
    uvec beta_order = stable_sort_index(beta2, "descend");
    beta2 = (beta2(beta_order)).eval();

    vec s(p);
    vec w(p);
    vec betax(p);

    uvec idx_i(p);
    uvec idx_j(p);

    uword k = 0;

    for (uword i = 0; i < p; i++) {
      idx_i(k) = i;
      idx_j(k) = i;
      s(k)     = beta2(i) - sigma*lambda(i)*step_size;
      w(k)     = s(k);

      while ((k > 0) && (w[k - 1] <= w(k))) {
        k--;
        idx_j(k)  = i;
        s(k)     += s(k + 1);
        w(k)      = s(k) / (i - idx_i(k) + 1.0);
      }
      k++;
    }

    for (uword j = 0; j < k; j++) {
      double d = std::max(w(j), 0.0);
      for (uword i = idx_i(j); i <= idx_j(j); i++) {
        betax(i) = d;
      }
    }

    // reset order
    betax(beta_order) = betax;

    // reset sign and return
    return (betax % beta_sign).eval();
  };

  double
  primal(const arma::mat& beta)
  {
    using namespace arma;
    return dot(sigma*lambda, sort(abs(beta), "descending"));
  }

  double
  infeasibility(const arma::mat& grad)
  {
    using namespace arma;

    vec grad_sorted = sort(abs(grad), "descending");
    return std::max(cumsum(grad_sorted - sigma*lambda).max(), 0.0);
  }

  void
  step(const arma::uword i)
  {
    // lambda paths are currently not implemented for SLOPE
  }

  arma::uword
  pathLength()
  {
    return 1;
  }

  double
  lambdaInfeas()
  {
    return lambda(0)*sigma;
  }
};

class Lasso : public Penalty {
public:
  double lambda;
  const arma::vec lambda_path;

  Lasso(const arma::vec& lambda_path)
        : lambda_path(lambda_path), lambda(lambda_path(0)) {}

  arma::mat
  eval(const arma::mat& beta, const double step_size)
  {
    using namespace arma;

    return sign(beta) % clamp(abs(beta) - step_size*lambda, 0.0, datum::inf);
  };

  double
  primal(const arma::mat& beta)
  {
    return lambda*arma::norm(beta, 1);
  }

  double
  infeasibility(const arma::mat& grad)
  {
    using namespace arma;

    vec grad_sorted = sort(abs(grad), "descending");
    return std::max(cumsum(grad_sorted - lambda).max(), 0.0);
  }

  void
  step(const arma::uword i)
  {
    if (i < lambda_path.n_elem)
      lambda = lambda_path(i);
  }

  arma::uword
  pathLength()
  {
    return lambda_path.n_elem;
  }

  double
  lambdaInfeas()
  {
    return lambda;
  }
};

// helper to choose penalty
inline
std::unique_ptr<Penalty>
setupPenalty(const Rcpp::S4& args)
{
  std::string name = args.slot("name");
  arma::vec lambda = args.slot("lambda");

  if (name == "slope") {
    double sigma = args.slot("sigma");
    return std::unique_ptr<SLOPE>(new SLOPE{sigma, lambda});
  }

  // else lasso
  double lambda_scale = args.slot("lambda_scale");

  lambda /= lambda_scale;

  return std::unique_ptr<Lasso>(new Lasso{lambda});
}
