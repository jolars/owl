#ifndef GOLEM_UTILS
#define GOLEM_UTILS

#include <RcppArmadillo.h>

inline
arma::vec
sigmoid(const arma::vec& x)
{
  return 1.0/(1.0 + arma::exp(-x));
}

inline
double
sigmoid(const double x)
{
  return 1.0/(1.0 + std::exp(-x));
}

//' Clamp a value to [min, max]
//' @param x value to clamp
//' @param min min
//' @param max max
//' @noRd
template <typename T>
inline
T
clamp(const T& x, const T& min, const T& max)
{
  return x > max ? max : (x < min ? min : x);
}

class ConvergenceCheck {
public:
  ConvergenceCheck(const arma::rowvec& intercept_old,
                   const arma::vec& beta_old,
                   const double tol)
                   : intercept_old(intercept_old),
                     beta_old(beta_old),
                     tol(tol) {}

  bool
  operator()(const arma::rowvec& intercept_new,
             const arma::vec& beta_new)
  {
    double max_change = std::max(abs(intercept_new - intercept_old).max(),
                                 abs(beta_new - beta_old).max());
    double max_size   = std::max(abs(intercept_new).max(),
                                 abs(beta_new).max());

    bool all_zero  = (max_size == 0.0) && (max_change == 0.0);
    bool no_change = (max_size != 0.0) && (max_change/max_size <= tol);

    beta_old = beta_new;

    return all_zero || no_change;
  }

private:
  arma::rowvec intercept_old;
  arma::vec beta_old;
  const double tol;
};

template <typename T>
void
preprocessFeatures(T& x,
                   arma::rowvec& x_center,
                   arma::rowvec& x_scale,
                   const std::string& standardize)
{
  using namespace arma;

  auto p = x.n_cols;

  x_center.set_size(p);
  x_scale.set_size(p);

  // TODO(jolars): adapt for sparse input later

  if (standardize == "both" || standardize == "features") {
    x_center = mean(x);

    for (decltype(p) j = 0; j < p; ++j)
      x_scale(j) = norm(x.col(j));

    x_scale.replace(0, 1);

    for (decltype(p) j = 0; j < p; ++j) {
      x.col(j) -= x_center(j);
      x.col(j) /= x_scale(j);
    }

  } else {
    x_center.zeros();
    x_scale.ones();
  }
}

inline
void
unstandardize(arma::cube& intercepts,
              arma::cube& betas,
              const arma::rowvec& x_center,
              const arma::rowvec& x_scale,
              const arma::rowvec& y_center,
              const arma::rowvec& y_scale,
              const bool fit_intercept)
{
  using namespace arma;

  uword p = betas.n_rows;
  uword m = betas.n_cols;

  for (uword k = 0; k < m; ++k) {
    double x_bar_beta_sum = 0.0;

    for (uword j = 0; j < p; ++j) {
      betas.tube(j, k) *= y_scale(k)/x_scale(j);
      x_bar_beta_sum += accu(x_center(j) * betas.tube(j, k));
    }

    if (fit_intercept)
      intercepts.tube(0, k) = intercepts.tube(0, k) * y_scale(k) + y_center(k) - x_bar_beta_sum;
  }
}

inline
arma::vec
logSeq(const double from,
       const double to,
       const unsigned n)
{
  return arma::exp(arma::linspace(std::log(from), std::log(to), n));
}

#endif /* GOLEM_UTILS */
