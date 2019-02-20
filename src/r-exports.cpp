#include <RcppArmadillo.h>
#include "penalties.h"

// [[Rcpp::export]]
arma::vec
prox_slope_cpp(const arma::vec& y, const Rcpp::List& args)
{
  SLOPE penalty{args};

  return penalty.eval(y, 1.0);
}
