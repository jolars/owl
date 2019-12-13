#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "families.h"
#include "utils.h"
#include "proxes.h"

using namespace arma;

class Penalty {
public:
  virtual
  mat
  eval(const mat& beta, const vec& lambda, const double shrinkage) = 0;

  virtual
  double
  primal(const mat& beta, const vec& lambda) = 0;

  virtual
  double
  infeasibility(const mat& gradient, const vec& lambda) = 0;

  virtual
  uvec
  activeSet(const std::unique_ptr<Family>& family,
            const mat& y,
            const mat& gradient_prev,
            const mat& pseudo_gradient_prev,
            const vec& norms,
            const vec& lambda,
            const vec& lambda_prev,
            const std::string screening_rule) = 0;

  // returns the indices of coefficients for which the kkt test fails
  uvec
  kktCheck(const mat&   gradient,
           const mat&   beta,
           const vec&   lambda,
           const double tol = 1e-6)
  {
    uvec nonzeros = find(beta != 0);
    uvec ord = sort_index(abs(gradient), "descend");
    vec abs_gradient_sorted = abs(gradient(ord));

    double rh = std::max(std::sqrt(datum::eps), tol*lambda(0));

    uvec out = cumsum(abs_gradient_sorted - lambda) <= rh;
    out(ord) = out;
    out(nonzeros).ones();

    return find(out == 0);
  }
};

class SLOPE : public Penalty {
public:
  mat
  eval(const mat& beta, const vec& lambda, const double shrinkage)
  {
    mat out = slopeProx(vectorise(beta), lambda, shrinkage);

    return reshape(out, size(beta));
  };

  double
  primal(const mat& beta, const vec& lambda)
  {
    return dot(lambda, sort(abs(vectorise(beta)), "descending"));
  }

  double
  infeasibility(const mat& gradient, const vec& lambda)
  {
    vec abs_gradient_sorted = sort(abs(vectorise(gradient)), "descending");
    return std::max(cumsum(abs_gradient_sorted - lambda).max(), 0.0);
  }

  uvec
  activeSet(const std::unique_ptr<Family>& family,
            const mat& y,
            const mat& gradient_prev,
            const mat& pseudo_gradient_prev,
            const vec& norms,
            const vec& lambda,
            const vec& lambda_prev,
            const std::string screening_rule);
};

class GroupSLOPE : public Penalty {
public:
  const field<uvec> group_id;
  const uword n_groups;

  GroupSLOPE(const field<uvec>& group_id)
             : group_id(group_id),
               n_groups(group_id.n_elem) {}

  mat
  eval(const mat& beta,
       const vec& lambda,
       const double shrinkage)
  {
    vec group_norms(n_groups);

    for (uword i = 0; i < n_groups; ++i)
      group_norms(i) = norm(beta(group_id(i)), "fro");

    auto prox_norms = slopeProx(group_norms, lambda, shrinkage);

    mat prox_solution(beta.n_rows*beta.n_cols, 1);

    for (uword i = 0; i < n_groups; ++i) {
      uvec idx = group_id(i);
      prox_solution(idx) = beta(idx) * (prox_norms(i)/group_norms(i));
    }

    return reshape(prox_solution, size(beta));
  };

  double
  primal(const mat& beta, const vec& lambda)
  {
    vec beta_norms(n_groups);

    for (uword i = 0; i < n_groups; ++i)
      beta_norms(i) = norm(beta(group_id(i)), "fro");

    return dot(lambda, sort(beta_norms, "descending"));
  }

  double
  infeasibility(const mat& gradient, const vec& lambda)
  {
    vec gradient_norms(n_groups);

    for (uword i = 0; i < n_groups; ++i)
      gradient_norms(i) = norm(gradient(group_id(i)), "fro");

    const vec gradient_norms_sorted = sort(gradient_norms, "descending");
    return std::max(cumsum(gradient_norms_sorted - lambda).max(), 0.0);
  }

  uvec
  activeSet(const std::unique_ptr<Family>& family,
            const mat& y,
            const mat& gradient_prev,
            const mat& pseudo_gradient_prev,
            const vec& norms,
            const vec& lambda,
            const vec& lambda_prev,
            const std::string screening_rule);
};

// helper to choose penalty
inline
std::unique_ptr<Penalty>
setupPenalty(const Rcpp::List& args, const Rcpp::List& groups)
{
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
