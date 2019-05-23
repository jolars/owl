#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "families.h"
#include "utils.h"
#include "proxes.h"

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
  const double sigma;
  const arma::vec lambda;

  SLOPE(const double sigma, const arma::vec& lambda)
        : sigma(sigma), lambda(lambda) {}

  arma::mat
  eval(const arma::mat& beta, const double step_size)
  {
    return slopeProx(beta, step_size, lambda, sigma);
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

class GroupSLOPE : public Penalty {
public:
  const double sigma;
  const arma::vec lambda;
  const arma::field<arma::uvec> groups;
  const arma::uword n_groups;

  GroupSLOPE(const double sigma,
             const arma::vec& lambda,
             const arma::field<arma::uvec>& groups)
             : sigma(sigma),
               lambda(lambda),
               groups(groups),
               n_groups(groups.n_elem) {}

  arma::mat
  eval(const arma::mat& beta, const double step_size)
  {
    using namespace arma;

    vec group_norms(n_groups);

    for (uword i = 0; i < n_groups; ++i)
      group_norms(i) = norm(beta(groups(i)), "fro");

    auto prox_norms = slopeProx(group_norms, step_size, lambda, sigma);

    vec prox_solution(size(beta));

    for (uword i = 0; i < n_groups; ++i) {
      uvec idx = groups(i);
      prox_solution(idx) = beta(idx) * (prox_norms(i)/group_norms(i));
    }

    return prox_solution;
  };

  double
  primal(const arma::mat& beta)
  {
    using namespace arma;

    vec beta_norms(n_groups);

    for (uword i = 0; i < n_groups; ++i)
      beta_norms(i) = norm(beta(groups(i)), "fro");

    return dot(sigma*lambda, sort(beta_norms, "descending"));
  }

  double
  infeasibility(const arma::mat& grad)
  {
    using namespace arma;

    vec grad_norms(n_groups);

    for (uword i = 0; i < n_groups; ++i)
      grad_norms(i) = norm(grad(groups(i)), "fro");

    const vec grad_norms_sorted = sort(grad_norms, "descending");
    return std::max(cumsum(grad_norms_sorted - sigma*lambda).max(), 0.0);
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
  const arma::vec lambda_path;
  double lambda;

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
    return std::max(cumsum(sort(abs(grad), "descending") - lambda).max(), 0.0);
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
  using namespace arma;

  std::string name = args.slot("name");
  vec lambda = args.slot("lambda");

  if (name == "slope") {

    double sigma = args.slot("sigma");
    return std::unique_ptr<SLOPE>(new SLOPE{sigma, lambda});

  } else if (name == "group_slope") {

    double sigma = args.slot("sigma");
    bool orthogonalize = args.slot("orthogonalize");

    field<uvec> groups = orthogonalize ? args.slot("ortho_group_id")
                                       : args.slot("group_id");

    groups.for_each([](uvec& x) {x -= 1;}); // fix indexing for c++

    return std::unique_ptr<GroupSLOPE>(new GroupSLOPE{sigma, lambda, groups});
  }

  // else lasso
  double lambda_scale = args.slot("lambda_scale");

  lambda /= lambda_scale;

  return std::unique_ptr<Lasso>(new Lasso{lambda});
}
