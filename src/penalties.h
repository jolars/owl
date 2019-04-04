#ifndef GOLEM_PENALTY_
#define GOLEM_PENALTY_

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
  dual(const arma::mat& grad, const arma::mat& beta) = 0;

  virtual
  double
  infeasibility(const arma::mat& grad) = 0;

  virtual
  void
  step() = 0;

  virtual
  Rcpp::List
  getParams(const arma::vec& y_scale) = 0;

  virtual
  arma::uword
  pathLength() = 0;
};

class SLOPE : public Penalty {
private:
  arma::vec lambda;
  double sigma;

public:
  SLOPE(const Rcpp::List& args)
  {
    using namespace arma;
    using Rcpp::as;

    lambda = as<arma::vec>(args["lambda"]);
    sigma = as<double>(args["sigma"]);

    auto lambda_type = as<std::string>(args["lambda_type"]);
    auto sigma_type  = as<std::string>(args["sigma_type"]);
    auto fdr         = as<double>(args["fdr"]);
    auto n           = as<uword>(args["n"]);
    auto p           = as<uword>(args["p"]);

    if (lambda_type == "bhq" || lambda_type == "gaussian") {
      // begin creating the BHQ sequence
      lambda.set_size(p);
      vec q = regspace<vec>(1, p) * fdr/(2.0*p);

      for (uword i = 0; i < p; ++i)
        lambda(i) = R::qnorm(1.0 - q(i), 0.0, 1.0, 1, 0);

      if (lambda_type == "gaussian") {
        // first entry corresponds to the BHQ sequence

        if (p >= 2) {
          double sum_sq = 0;
          for (uword i = 1; i < p; ++i) {
            sum_sq += std::pow(lambda(i - 1), 2);
            double w = std::max(1.0, static_cast<double>(n - i - 1));
            lambda(i) *= std::sqrt(1.0 + sum_sq/w);
          }
        }

        auto lambda_min_index = lambda.index_min();
        auto lambda_min = lambda.min();
        lambda.tail(p - lambda_min_index).fill(lambda_min);
      }
    }

    if (lambda.n_elem != p)
      Rcpp::stop("lambda sequence must be as long as there are variables");

  };

  arma::mat
  eval(const arma::mat& beta, const double L)
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
      s(k)     = beta2(i) - sigma*lambda(i)/L;
      w(k)     = s(k);

      while ((k > 0) && (w[k-1] <= w(k))) {
        k--;
        idx_j(k)  = i;
        s(k)     += s(k+1);
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
  dual(const arma::mat& grad, const arma::mat& beta)
  {
    using namespace arma;
    return dot(grad, log(grad));
  }

  double
  infeasibility(const arma::mat& grad)
  {
    using namespace arma;

    vec grad_sorted = sort(abs(grad), "descending");
    return std::max(cumsum(grad_sorted - sigma*lambda).max(), 0.0);
  }

  void
  step()
  {
    // does nothing for slope
  }

  arma::uword
  pathLength()
  {
    return 1;
  }

  Rcpp::List
  getParams(const arma::vec& y_scale)
  {
    using namespace Rcpp;

    return List::create(Named("name")   = "slope",
                        Named("lambda") = wrap(lambda),
                        Named("sigma")  = sigma);
  }

};

class ElasticNet : public Penalty {
private:
  arma::vec lambda;
  double alpha;
  arma::uword lambda_ind = 0;

public:
  ElasticNet(const Rcpp::List& args,
             const arma::mat& x,
             const arma::mat& y,
             const arma::rowvec& y_scale,
             const std::unique_ptr<Family>& family)
  {
    using namespace arma;
    using Rcpp::as;

    lambda                = as<vec>(args["lambda"]);
    alpha                 = as<double>(args["alpha"]);
    auto n_lambda         = as<uword>(args["n_lambda"]);
    auto lambda_min_ratio = as<double>(args["lambda_min_ratio"]);
    auto lambda_type      = as<std::string>(args["lambda_type"]);

    if (lambda_type == "auto") {
      double lambda_max =
        family->lambdaMax(x, y, y_scale)/std::max(alpha, 0.001);

      lambda = logSeq(lambda_max, lambda_max*lambda_min_ratio, n_lambda);
    }

    lambda /= y_scale.max();
  };

  arma::mat
  eval(const arma::mat& beta, const double L)
  {
    using namespace arma;

    return sign(beta)
           % clamp(abs(beta) - alpha*lambda[lambda_ind]*(1.0/L), 0.0, datum::inf);
  };

  double
  primal(const arma::mat& beta)
  {
    using namespace arma;
    return alpha*lambda[lambda_ind]*accu(abs(beta));
  }

  double
  dual(const arma::mat& grad, const arma::mat& beta)
  {
    return arma::as_scalar(grad.t() * arma::log(grad));
  }

  double
  infeasibility(const arma::mat& grad)
  {
    using namespace arma;

    vec grad_sorted = sort(abs(grad), "descending");
    return std::max(cumsum(grad_sorted - alpha*lambda[lambda_ind]).max(), 0.0);
  }

  void
  step()
  {
    lambda_ind++;
  }

  arma::uword
  pathLength()
  {
    return lambda.n_elem;
  }

  Rcpp::List
  getParams(const arma::vec& y_scale)
  {
    using namespace Rcpp;

    arma::vec lambda_out = lambda*y_scale.max();

    return List::create(Named("name")   = "elasticNet",
                        Named("lambda") = wrap(lambda_out),
                        Named("alpha")  = alpha);
  }
};

// helper to choose penalty
inline
std::unique_ptr<Penalty>
setupPenalty(const Rcpp::List& args,
             const arma::mat& x,
             const arma::mat& y,
             const arma::rowvec& y_scale,
             const std::unique_ptr<Family>& family)
{
  auto name = Rcpp::as<std::string>(args["name"]);

  if (name == "slope") {
    return std::unique_ptr<SLOPE>(new SLOPE{args});
  }

  return
    std::unique_ptr<ElasticNet>(new ElasticNet{args, x, y, y_scale, family});
}


#endif /* GOLEM_PENALTIES_ */
