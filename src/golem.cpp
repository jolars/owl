#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List
denseGolem(arma::mat X,
           arma::vec y,
           Rcpp::List control)
{
  return control;
}
