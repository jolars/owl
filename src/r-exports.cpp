#include <RcppArmadillo.h>
#include "penalties.h"
#include "utils.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec
prox_slope_cpp(const arma::mat& y, const Rcpp::List& args)
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
lambdaMax(const T& x,
          const arma::mat& y,
          const arma::vec& x_center,
          const arma::vec& x_scale,
          const arma::vec& y_scale,
          const uword n_targets,
          const std::string& family,
          const bool standardize_features,
          const bool is_sparse)
{
  const uword p = x_center.n_elem;
  mat lambda_max(p, n_targets);

  if (family == "binomial") {
    vec y_new = (y + 1)/2;

    // standardize
    double y_center = mean(y_new);
    y_new -= y_center;

    lambda_max = x.t() * y_new;

    if (is_sparse && standardize_features) {
      for (uword j = 0; j < p; ++j)
        lambda_max(j) -= accu(y_new * x_center(j)/x_scale(j));
    }

  } else if (family == "multinomial") {

    uword n = y.n_rows;

    uvec y_classes = conv_to<uvec>::from(y + 0.1);
    mat y_map(n, n_targets);

    // for (uword i = 0; i < n; ++i) {
    //   auto c = static_cast<uword>(y(i) + 0.1);
    //   y_map(i, c) = 1.0;
    // }

    for (uword k = 0; k < n_targets; ++k) {
      y_map.col(k) = conv_to<vec>::from(y_classes == k);
    }

    rowvec y_bar = mean(y_map);
    rowvec y_std = stddev(y_map, 1);

    for (uword k = 0; k < n_targets; ++k) {
      y_map.col(k) -= y_bar(k);
      y_map.col(k) /= y_std(k);
    }

    lambda_max = x.t() * y_map;

    for (uword k = 0; k < n_targets; ++k) {
      lambda_max.col(k) *= y_std(k);
    }

  } else {

    lambda_max = x.t() * y;

    if (is_sparse && standardize_features) {
      for (uword j = 0; j < p; ++j)
        lambda_max(j) -= accu(y * x_center(j)/x_scale(j));
    }
  }

  return abs(vectorise(lambda_max));
}

// [[Rcpp::export]]
arma::vec
lambdaMax(SEXP x,
          const arma::mat& y,
          const arma::vec& x_center,
          const arma::vec& x_scale,
          const arma::vec& y_scale,
          const arma::uword n_targets,
          const std::string& family,
          const bool standardize_features)
{
  vec lambda_max;
  const bool is_sparse = isSparse(x);

  if (is_sparse)
    return lambdaMax(Rcpp::as<sp_mat>(x),
                     y,
                     x_center,
                     x_scale,
                     y_scale,
                     n_targets,
                     family,
                     standardize_features,
                     is_sparse);
  else
    return lambdaMax(Rcpp::as<mat>(x),
                     y,
                     x_center,
                     x_scale,
                     y_scale,
                     n_targets,
                     family,
                     standardize_features,
                     is_sparse);
}

// template <typename T>
// arma::vec
// rowNorms(const T& x, const arma::uword norm_type = 2)
// {
//   using namespace arma;
//
//   uword p = x.n_cols;
//
//   vec norms(p);
//
//   for (decltype(p) i = 0; i < p; ++i)
//     norms(i) = norm(x.row(i), norm_type);
//
//   return norms;
// }

// // [[Rcpp::export]]
// arma::vec
// rowNormsSparse(const arma::sp_mat& x,
//                const arma::uword norm_type = 2)
// {
//   return rowNorms(x, norm_type);
// }
//
// // [[Rcpp::export]]
// arma::vec
// rowNormsDense(const arma::mat& x, const arma::uword norm_type = 2)
// {
//   return rowNorms(x, norm_type);
// }


// double
// maxSquaredRowNorm(const arma::mat& x,
//                   const arma::rowvec& x_scaled_center,
//                   const bool standardize_features)
// {
//   return arma::norm(arma::square(x), "inf");
// }
//
// double
// maxSquaredRowNorm(const arma::sp_mat& x,
//                   const arma::rowvec& x_scaled_center,
//                   const bool standardize_features)
// {
//   using namespace arma;
//
//   double max_squared_row_norms = 0.0;
//
//   if (standardize_features) {
//     const uword n = x.n_rows;
//
//     vec squared_row_norms(n);
//
//     for (uword i = 0; i < n; ++i)
//       squared_row_norms(i) = accu(square(x.row(i) - x_scaled_center));
//
//     max_squared_row_norms = squared_row_norms.max();
//
//   } else {
//     max_squared_row_norms = sum(square(x), 1).max();
//   }
//
//   return max_squared_row_norms;
// }
//
// // [[Rcpp::export]]
// double
// maxSquaredRowNorm(SEXP x,
//                   const arma::rowvec& x_scaled_center,
//                   const bool standardize_features)
// {
//   if (Rf_isS4(x)) {
//     if (Rf_inherits(x, "dgCMatrix"))
//       return maxSquaredRowNorm(Rcpp::as<arma::sp_mat>(x),
//                                x_scaled_center,
//                                standardize_features);
//   }
//
//   return maxSquaredRowNorm(Rcpp::as<arma::mat>(x),
//                            x_scaled_center,
//                            standardize_features);
// }

Rcpp::List
standardizeDense(arma::mat x)
{
  const uword p = x.n_cols;
  const uword n = x.n_rows;

  vec x_center(p);
  vec x_scale(p);

  for (uword j = 0; j < p; ++j) {
    x_center(j) = accu(x.col(j))/n;
    x_scale(j) = norm(x.col(j) - x_center(j), 2);
  }

  // don't scale zero-variance predictors
  x_scale.replace(0, 1);

  for (uword j = 0; j < p; ++j) {
    x.col(j) -= x_center(j);
    x.col(j) /= x_scale(j);
  }

  return List::create(
    Named("x") = wrap(x),
    Named("x_center") = wrap(x_center),
    Named("x_scale") = wrap(x_scale)
  );
}

Rcpp::List
standardizeSparse(arma::sp_mat x)
{
  const uword p = x.n_cols;
  const uword n = x.n_rows;

  vec x_center(p);
  vec x_scale(p);

  for (uword j = 0; j < p; ++j) {
    x_center(j) = accu(x.col(j))/n;
    x_scale(j) = norm(x.col(j) - x_center(j), 2);
  }

  // don't scale zero-variance predictors
  x_scale.replace(0, 1);

  for (uword j = 0; j < p; ++j) {
    x.col(j) /= x_scale(j);
  }

  return List::create(
    Named("x") = wrap(x),
    Named("x_center") = wrap(x_center),
    Named("x_scale") = wrap(x_scale)
  );
}

// [[Rcpp::export]]
Rcpp::List
standardize(SEXP x)
{
  bool is_sparse = isSparse(x);

  if (is_sparse)
    return standardizeSparse(Rcpp::as<sp_mat>(x));
  else
    return standardizeDense(Rcpp::as<mat>(x));
}

// [[Rcpp::export]]
arma::mat
tester(arma::uvec y, arma::mat lin_pred)
{

  mat bla = lin_pred.each_row([](const rowvec& a) {return log(accu(exp(a - a.max()))) + a.max();});

  return bla;
}
