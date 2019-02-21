#ifndef GOLEM_PENALTY_
#define GOLEM_PENALTY_

#include <RcppArmadillo.h>
#include <memory>

class Penalty {
public:
  virtual
  arma::vec
  eval(arma::vec y, const double L) = 0;

  virtual
  double
  primal(const arma::vec& beta) = 0;

  virtual
  double
  dual(const arma::vec& grad, const arma::vec& beta) = 0;

  virtual
  double
  infeasibility(const arma::vec& grad) = 0;

  virtual
  double
  getSigma() = 0;

  virtual
  void
  setSigma(const double) = 0;

  virtual
  arma::vec
  getLambda() = 0;
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

  arma::vec
  eval(arma::vec beta, const double L)
  {
    using namespace arma;

    uword p = beta.n_elem;

    // collect sign of beta and work with sorted absolutes
    vec beta_sign = sign(beta);
    beta = abs(beta);
    uvec beta_order = stable_sort_index(beta, "descend");
    beta = (beta(beta_order)).eval();

    vec s(p);
    vec w(p);
    vec betax(p);

    uvec idx_i(p);
    uvec idx_j(p);

    uword k = 0;

    for (uword i = 0; i < p; i++) {
      idx_i(k) = i;
      idx_j(k) = i;
      s(k)     = beta(i) - sigma*lambda(i)/L;
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
  primal(const arma::vec& beta)
  {
    using namespace arma;
    return dot(sigma*lambda, sort(abs(beta), "descending"));
  }

  double
  dual(const arma::vec& grad, const arma::vec& beta)
  {
    using namespace arma;
    return dot(grad, log(grad));
  }

  double
  infeasibility(const arma::vec& grad)
  {
    using namespace arma;

    vec grad_sorted = sort(abs(grad), "descending");
    return std::max(cumsum(grad_sorted - sigma*lambda).max(), 0.0);
  }

  double
  getSigma()
  {
    return sigma;
  }

  void
  setSigma(const double sigma_new)
  {
    sigma = sigma_new;
  }

  arma::vec
  getLambda()
  {
    return lambda;
  }
};


// helper to choose penalty
inline
std::unique_ptr<Penalty>
setupPenalty(const Rcpp::List& args)
{
  //auto name = Rcpp::as<std::string>(args["name"]);

  return std::unique_ptr<SLOPE>(new SLOPE{args});
}


#endif /* GOLEM_PENALTIES_ */
