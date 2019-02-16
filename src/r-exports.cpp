#include <RcppArmadillo.h>
#include "penalties.h"

// [[Rcpp::export]]
arma::vec
prox_slope_cpp(arma::vec y,
               const arma::vec lambda)
{
  SLOPE penalty(lambda);

  return penalty.eval(y, 1.0);
}
