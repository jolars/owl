#include <RcppArmadillo.h>
#include "proxes.h"

// [[Rcpp::export]]
arma::vec
prox_slope_cpp(arma::vec y,
               const arma::vec lambda)
{
  SLOPE prox(lambda);

  return prox.eval(y, 1.0);
}
