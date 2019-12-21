#include <RcppArmadillo.h>
#include "prox.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec prox_slope_cpp(const arma::vec& y, const Rcpp::List& args)
{
  auto sigma = Rcpp::as<arma::vec>(args["sigma"]);
  auto lambda = Rcpp::as<arma::vec>(args["lambda"]);

  return prox(y, lambda*sigma(0));
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

