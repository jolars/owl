#include <RcppArmadillo.h>
#include "penalties.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec prox_slope_cpp(const arma::mat& y, const Rcpp::List& args)
{
  auto sigma = Rcpp::as<arma::vec>(args["sigma"]);
  auto lambda = Rcpp::as<arma::vec>(args["lambda"]);

  SLOPE penalty;

  return penalty.eval(y, lambda*sigma(0), 1.0);
}

//
// // [[Rcpp::export]]
// arma::mat
// tester(arma::uvec y, arma::mat lin_pred)
// {
//
//   mat bla = lin_pred.each_row([](const rowvec& a) {return log(accu(exp(a - a.max()))) + a.max();});
//
//   return bla;
// }

