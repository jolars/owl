#pragma once

#include <RcppArmadillo.h>
#include <memory>
#include "../penalties.h"
#include "../families.h"
#include "../utils.h"

using namespace arma;

struct Results {
  rowvec intercept;
  mat beta;
  uword passes;
  std::vector<double> primals;
  std::vector<double> duals;
  std::vector<double> infeasibilities;
  std::vector<double> time;
  std::vector<unsigned> line_searches;

  Results() {}

  Results(rowvec intercept,
          mat beta,
          uword passes,
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
  fit(const sp_mat& x,
      const mat& y,
      const std::unique_ptr<Family>& family,
      const std::unique_ptr<Penalty>& penalty,
      const rowvec& intercept_init,
      const mat& beta_init,
      const bool fit_intercept,
      const vec& lambda,
      const vec& x_center,
      const vec& x_scale) = 0;

  virtual
  Results
  fit(const mat& x,
      const mat& y,
      const std::unique_ptr<Family>& family,
      const std::unique_ptr<Penalty>& penalty,
      const rowvec& intercept_init,
      const mat& beta_init,
      const bool fit_intercept,
      const vec& lambda,
      const vec& x_center,
      const vec& x_scale) = 0;
};

