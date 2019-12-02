#include <RcppArmadillo.h>
#include "penalties.h"
#include "utils.h"

// [[Rcpp::export]]
arma::vec
prox_slope_cpp(const arma::vec& y, const Rcpp::List& args)
{
  auto sigma = Rcpp::as<arma::vec>(args["sigma"]);
  auto lambda = Rcpp::as<arma::vec>(args["lambda"]);

  SLOPE penalty;

  return penalty.eval(y, lambda*sigma(0), 1.0);
}


// [[Rcpp::export]]
arma::vec
standardizedSparseColNorms(const arma::sp_mat& x,
                           const arma::vec& x_center)
{
  using namespace arma;

  const uword p = x.n_cols;
  vec norms(p, fill::zeros);

  for (uword j = 0; j < p; ++j)
    norms(j) = norm(x.col(j) - x_center(j), 2);

  return norms;
}

template <typename T>
arma::vec
colNorms(const T& x, const arma::uword norm_type = 2)
{
  using namespace arma;

  uword p = x.n_cols;

  arma::vec norms(p);

  for (decltype(p) i = 0; i < p; ++i)
    norms(i) = arma::norm(x.col(i), norm_type);

  return norms;
}

// [[Rcpp::export]]
arma::vec
colNormsSparse(const arma::sp_mat& x, const arma::uword norm_type = 2)
{
  return colNorms(x, norm_type);
}

// [[Rcpp::export]]
arma::vec
colNormsDense(const arma::mat& x, const arma::uword norm_type = 2)
{
  return colNorms(x, norm_type);
}

template <typename T>
arma::vec
rowNorms(const T& x, const arma::uword norm_type = 2)
{
  using namespace arma;

  uword p = x.n_cols;

  vec norms(p);

  for (decltype(p) i = 0; i < p; ++i)
    norms(i) = norm(x.row(i), norm_type);

  return norms;
}

// [[Rcpp::export]]
arma::vec
rowNormsSparse(const arma::sp_mat& x,
               const arma::uword norm_type = 2)
{
  return rowNorms(x, norm_type);
}

// [[Rcpp::export]]
arma::vec
rowNormsDense(const arma::mat& x, const arma::uword norm_type = 2)
{
  return rowNorms(x, norm_type);
}


double
maxSquaredRowNorm(const arma::mat& x,
                  const arma::rowvec& x_scaled_center,
                  const bool standardize_features)
{
  return arma::norm(arma::square(x), "inf");
}

double
maxSquaredRowNorm(const arma::sp_mat& x,
                  const arma::rowvec& x_scaled_center,
                  const bool standardize_features)
{
  using namespace arma;

  double max_squared_row_norms = 0.0;

  if (standardize_features) {
    const uword n = x.n_rows;

    vec squared_row_norms(n);

    for (uword i = 0; i < n; ++i)
      squared_row_norms(i) = accu(square(x.row(i) - x_scaled_center));

    max_squared_row_norms = squared_row_norms.max();

  } else {
    max_squared_row_norms = sum(square(x), 1).max();
  }

  return max_squared_row_norms;
}

// [[Rcpp::export]]
double
maxSquaredRowNorm(SEXP x,
                  const arma::rowvec& x_scaled_center,
                  const bool standardize_features)
{
  if (Rf_isS4(x)) {
    if (Rf_inherits(x, "dgCMatrix"))
      return maxSquaredRowNorm(Rcpp::as<arma::sp_mat>(x),
                               x_scaled_center,
                               standardize_features);
  }

  return maxSquaredRowNorm(Rcpp::as<arma::mat>(x),
                           x_scaled_center,
                           standardize_features);
}




