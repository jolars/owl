#pragma once

#include <RcppArmadillo.h>

using namespace Rcpp;
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
  double deviance;

  Results() {}

  Results(rowvec intercept,
          mat beta,
          uword passes,
          std::vector<double> primals,
          std::vector<double> duals,
          std::vector<double> infeasibilities,
          std::vector<double> time,
          std::vector<unsigned> line_searches,
          double deviance)
    : intercept(intercept),
      beta(beta),
      passes(passes),
      primals(primals),
      duals(duals),
      infeasibilities(infeasibilities),
      time(time),
      line_searches(line_searches),
      deviance(deviance) {}
};
