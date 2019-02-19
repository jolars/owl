#include <RcppArmadillo.h>
#include "penalties.h"

// [[Rcpp::export]]
arma::vec
prox_slope_cpp(arma::vec y,
               const arma::vec lambda)
{
  SLOPE penalty(lambda, "user", 1, "user", y.n_elem, 1, 0.2);

  return penalty.eval(y, 1.0);
}
