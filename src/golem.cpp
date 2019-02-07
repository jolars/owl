#include <RcppArmadillo.h>
using namespace Rcpp;

//' Fit a Model with sgdnet
//'
//' This main use of this function is calling the templated SetupSgdnet()
//' so that the dense and sparse implementations are compiled and
//' called appropriately. The control parameters in `control` are just
//' passed along.
//'
//' @param x a double
//' @param y a double
//'
//' @return Something
// [[Rcpp::export]]
double what(double x, double y)
{
  return x*y;
}
