#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "families.h"
#include "utils.h"
#include "proxes.h"

class Penalty {
public:
  virtual
  arma::vec
  eval(const arma::vec& beta,
       const arma::vec& lambda,
       const double shrinkage) = 0;

  virtual
  double
  primal(const arma::vec& beta, const arma::vec& lambda) = 0;

  virtual
  double
  infeasibility(const arma::vec& gradient, const arma::vec& lambda) = 0;

  virtual
  arma::uvec
  activeSet(const std::unique_ptr<Family>& family,
            const arma::vec& y,
            const arma::vec& gradient_prev,
            const arma::vec& pseudo_gradient_prev,
            const arma::vec& norms,
            const arma::vec& lambda,
            const arma::vec& lambda_prev,
            const std::string screening_rule) = 0;
};

class SLOPE : public Penalty {
public:
  arma::vec
  eval(const arma::vec& beta,
       const arma::vec& lambda,
       const double shrinkage)
  {
    return slopeProx(beta, lambda, shrinkage);
  };

  double
  primal(const arma::vec& beta, const arma::vec& lambda)
  {
    using namespace arma;
    return dot(lambda, sort(abs(beta), "descending"));
  }

  double
  infeasibility(const arma::vec& gradient, const arma::vec& lambda)
  {
    using namespace arma;

    vec gradient_sorted = sort(abs(gradient), "descending");
    return std::max(cumsum(gradient_sorted - lambda).max(), 0.0);
  }

  arma::uvec
  activeSet(const std::unique_ptr<Family>& family,
            const arma::vec& y,
            const arma::vec& gradient_prev,
            const arma::vec& pseudo_gradient_prev,
            const arma::vec& norms,
            const arma::vec& lambda,
            const arma::vec& lambda_prev,
            const std::string screening_rule);
};

class GroupSLOPE : public Penalty {
public:
  const arma::field<arma::uvec> group_id;
  const arma::uword n_groups;

  GroupSLOPE(const arma::field<arma::uvec>& group_id)
             : group_id(group_id),
               n_groups(group_id.n_elem) {}

  arma::vec
  eval(const arma::vec& beta,
       const arma::vec& lambda,
       const double shrinkage)
  {
    using namespace arma;

    vec group_norms(n_groups);

    for (uword i = 0; i < n_groups; ++i)
      group_norms(i) = norm(beta(group_id(i)), "fro");

    auto prox_norms = slopeProx(group_norms, lambda, shrinkage);

    vec prox_solution(size(beta));

    for (uword i = 0; i < n_groups; ++i) {
      uvec idx = group_id(i);
      prox_solution(idx) = beta(idx) * (prox_norms(i)/group_norms(i));
    }

    return prox_solution;
  };

  double
  primal(const arma::vec& beta, const arma::vec& lambda)
  {
    using namespace arma;

    vec beta_norms(n_groups);

    for (uword i = 0; i < n_groups; ++i)
      beta_norms(i) = norm(beta(group_id(i)), "fro");

    return dot(lambda, sort(beta_norms, "descending"));
  }

  double
  infeasibility(const arma::vec& gradient, const arma::vec& lambda)
  {
    using namespace arma;

    vec gradient_norms(n_groups);

    for (uword i = 0; i < n_groups; ++i)
      gradient_norms(i) = norm(gradient(group_id(i)), "fro");

    const vec gradient_norms_sorted = sort(gradient_norms, "descending");
    return std::max(cumsum(gradient_norms_sorted - lambda).max(), 0.0);
  }

  arma::uvec
  activeSet(const std::unique_ptr<Family>& family,
            const arma::vec& y,
            const arma::vec& gradient_prev,
            const arma::vec& pseudo_gradient_prev,
            const arma::vec& norms,
            const arma::vec& lambda,
            const arma::vec& lambda_prev,
            const std::string screening_rule);
};

// helper to choose penalty
inline
std::unique_ptr<Penalty>
setupPenalty(const Rcpp::List& args, const Rcpp::List& groups)
{
  using namespace arma;
  using Rcpp::as;

  auto name = as<std::string>(args["name"]);

  if (name == "group_slope") {

    auto sigma = as<vec>(args["sigma"]);
    bool orthogonalize = as<bool>(groups["orthogonalize"]);

    field<uvec> group_id =
      orthogonalize ? as<field<uvec>>(groups["ortho_group_id"])
                    : as<field<uvec>>(groups["group_id"]);

    group_id.for_each([](uvec& x) {x -= 1;}); // fix indexing for c++

    return std::unique_ptr<GroupSLOPE>(new GroupSLOPE{group_id});
  }

  return std::unique_ptr<SLOPE>(new SLOPE);
}
