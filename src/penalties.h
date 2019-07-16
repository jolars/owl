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
  eval(const arma::mat& y,
       const double L,
       const arma::uword k) = 0;

  virtual
  double
  primal(const arma::mat& beta, const arma::uword k) = 0;

  virtual
  double
  infeasibility(const arma::mat& grad, const arma::uword k) = 0;

  virtual
  double
  lambdaInfeas(const arma::uword k) = 0;
};

class SLOPE : public Penalty {
public:
  const arma::vec sigma;
  const arma::vec lambda;

  SLOPE(const arma::vec& sigma, const arma::vec& lambda)
        : sigma(sigma), lambda(lambda) {}

  arma::mat
  eval(const arma::mat& beta,
       const double step_size,
       const arma::uword k)
  {
    return slopeProx(beta, step_size, lambda, sigma(k));
  };

  double
  primal(const arma::mat& beta, const arma::uword k)
  {
    using namespace arma;
    return dot(sigma(k)*lambda, sort(abs(beta), "descending"));
  }

  double
  infeasibility(const arma::mat& grad, const arma::uword k)
  {
    using namespace arma;

    vec grad_sorted = sort(abs(grad), "descending");
    return std::max(cumsum(grad_sorted - sigma(k)*lambda).max(), 0.0);
  }

  double
  lambdaInfeas(const arma::uword k)
  {
    return lambda(0)*sigma(k);
  }
};

class GroupSLOPE : public Penalty {
public:
  const arma::vec sigma;
  const arma::vec lambda;
  const arma::field<arma::uvec> group_id;
  const arma::uword n_groups;

  GroupSLOPE(const arma::vec& sigma,
             const arma::vec& lambda,
             const arma::field<arma::uvec>& group_id)
             : sigma(sigma),
               lambda(lambda),
               group_id(group_id),
               n_groups(group_id.n_elem) {}

  arma::mat
  eval(const arma::mat& beta,
       const double step_size,
       const arma::uword k)
  {
    using namespace arma;

    vec group_norms(n_groups);

    for (uword i = 0; i < n_groups; ++i)
      group_norms(i) = norm(beta(group_id(i)), "fro");

    auto prox_norms =
      slopeProx(group_norms, step_size, lambda, sigma(k));

    vec prox_solution(size(beta));

    for (uword i = 0; i < n_groups; ++i) {
      uvec idx = group_id(i);
      prox_solution(idx) = beta(idx) * (prox_norms(i)/group_norms(i));
    }

    return prox_solution;
  };

  double
  primal(const arma::mat& beta, const arma::uword k)
  {
    using namespace arma;

    vec beta_norms(n_groups);

    for (uword i = 0; i < n_groups; ++i)
      beta_norms(i) = norm(beta(group_id(i)), "fro");

    return dot(sigma(k)*lambda, sort(beta_norms, "descending"));
  }

  double
  infeasibility(const arma::mat& grad, const arma::uword k)
  {
    using namespace arma;

    vec grad_norms(n_groups);

    for (uword i = 0; i < n_groups; ++i)
      grad_norms(i) = norm(grad(group_id(i)), "fro");

    const vec grad_norms_sorted = sort(grad_norms, "descending");
    return std::max(cumsum(grad_norms_sorted - sigma(k)*lambda).max(), 0.0);
  }

  double
  lambdaInfeas(const arma::uword k)
  {
    return lambda(0)*sigma(k);
  }
};

class Lasso : public Penalty {
public:
  const arma::vec lambda;

  Lasso(const arma::vec& lambda)
        : lambda(lambda) {}

  arma::mat
  eval(const arma::mat& beta,
       const double step_size,
       const arma::uword k)
  {
    using namespace arma;

    return sign(beta) % clamp(abs(beta) - step_size*lambda(k),
                0.0,
                datum::inf);
  };

  double
  primal(const arma::mat& beta, const arma::uword k)
  {
    return lambda(k)*arma::norm(beta, 1);
  }

  double
  infeasibility(const arma::mat& grad, const arma::uword k)
  {
    using namespace arma;
    return std::max(
      cumsum(sort(abs(grad), "descending") - lambda(k)).max(), 0.0);
  }

  double
  lambdaInfeas(const arma::uword k)
  {
    return lambda(k);
  }
};

// helper to choose penalty
inline
std::unique_ptr<Penalty>
setupPenalty(const Rcpp::List& args, const Rcpp::List& groups)
{
  using namespace arma;
  using Rcpp::as;

  auto name = as<std::string>(args["name"]);
  auto lambda = as<vec>(args["lambda"]);

  if (name == "slope") {

    auto sigma = as<vec>(args["sigma"]);
    return std::unique_ptr<SLOPE>(new SLOPE{sigma, lambda});

  } else if (name == "group_slope") {

    auto sigma = as<vec>(args["sigma"]);
    bool orthogonalize = as<bool>(groups["orthogonalize"]);

    field<uvec> group_id =
      orthogonalize ? as<field<uvec>>(groups["ortho_group_id"])
                    : as<field<uvec>>(groups["group_id"]);

    group_id.for_each([](uvec& x) {x -= 1;}); // fix indexing for c++

    return std::unique_ptr<GroupSLOPE>(new GroupSLOPE{sigma, lambda, group_id});
  }

  // else lasso
  double lambda_scale = as<double>(args["lambda_scale"]);

  lambda /= lambda_scale;

  return std::unique_ptr<Lasso>(new Lasso{lambda});
}
