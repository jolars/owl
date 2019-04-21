#include <RcppArmadillo.h>
#include "penalties.h"

// [[Rcpp::export]]
arma::vec
prox_slope_cpp(const arma::vec& y, const Rcpp::List& args)
{
  double sigma = Rcpp::as<double>(args["sigma"]);
  arma::vec lambda = Rcpp::as<arma::vec>(args["lambda"]);

  SLOPE penalty{sigma, lambda};

  return penalty.eval(y, 1.0);
}

// [[Rcpp::export]]
arma::vec
colNorms(const arma::mat& x, const arma::uword norm_type = 2)
{
  using namespace arma;

  uword p = x.n_cols;
  uword n = x.n_rows;

  vec norms(p);

  for (decltype(p) i = 0; i < p; ++i) {
    norms(i) = norm(x.col(i), norm_type);
  }

  return norms;
}

// [[Rcpp::export]]
arma::vec
rowNorms(const arma::mat& x, const arma::uword norm_type = 2)
{
  using namespace arma;

  uword p = x.n_cols;
  uword n = x.n_rows;

  vec norms(p);

  for (decltype(p) i = 0; i < p; ++i) {
    norms(i) = norm(x.row(i), norm_type);
  }

  return norms;
}
