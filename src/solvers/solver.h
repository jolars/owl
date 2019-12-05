#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "../penalties.h"
#include "../families.h"
#include "../utils.h"

struct Results {
  double intercept;
  arma::vec beta;
  arma::uword passes;
  std::vector<double> primals;
  std::vector<double> duals;
  std::vector<double> infeasibilities;
  std::vector<double> time;
  std::vector<unsigned> line_searches;

  Results() {}

  Results(double intercept,
          arma::vec beta,
          arma::uword passes,
          std::vector<double> primals,
          std::vector<double> duals,
          std::vector<double> infeasibilities,
          std::vector<double> time,
          std::vector<unsigned> line_searches)
          : intercept(intercept),
            beta(beta),
            passes(passes),
            primals(primals),
            duals(duals),
            infeasibilities(infeasibilities),
            time(time),
            line_searches(line_searches) {}
};

class Solver {
protected:
  bool diagnostics;

public:
  virtual
  Results
  fit(const arma::sp_mat& x,
      const arma::vec& y,
      const std::unique_ptr<Family>& family,
      const std::unique_ptr<Penalty>& penalty,
      const double intercept_init,
      const arma::vec& beta_init,
      const bool fit_intercept,
      const double lipschitz_constant,
      const arma::vec& lambda,
      const arma::vec& x_center,
      const arma::vec& x_scale) = 0;

  virtual
  Results
  fit(const arma::mat& x,
      const arma::vec& y,
      const std::unique_ptr<Family>& family,
      const std::unique_ptr<Penalty>& penalty,
      const double intercept_init,
      const arma::vec& beta_init,
      const bool fit_intercept,
      const double lipschitz_constant,
      const arma::vec& lambda,
      const arma::vec& x_center,
      const arma::vec& x_scale) = 0;
};

