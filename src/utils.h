#ifndef GOLEM_UTILS
#define GOLEM_UTILS

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
  ConvergenceCheck(const arma::vec& beta_old, const double tol)
    : beta_old(beta_old), tol(tol) {}

  bool
  operator()(const arma::vec& beta_new)
  {
    double max_change = (abs(beta_new - beta_old)).max();
    double max_size   = (abs(beta_new)).max();

    bool all_zero  = (max_size == 0.0) && (max_change == 0.0);
    bool no_change = (max_size != 0.0) && (max_change/max_size <= tol);

    beta_old = beta_new;

    return all_zero || no_change;
  }

private:
  arma::vec beta_old;
  const double tol;
};

template <typename T>
void
preprocessFeatures(T& X,
                   arma::rowvec& X_center,
                   arma::rowvec& X_scale,
                   const std::string& standardize)
{
  using namespace arma;

  auto p = X.n_cols;

  X_center.set_size(p);
  X_scale.set_size(p);

  // TODO(jolars): adapt for sparse input later

  if (standardize == "both" || standardize == "features") {
    X_center = mean(X);
    X_scale = stddev(X);
    X_scale.replace(0, 1); // avoid scaling by 0

    for (decltype(p) j = 0; j < p; ++j) {
      X.col(j) -= X_center(j);
      X.col(j) /= X_scale(j);
    }
  } else {
    X_center.zeros();
    X_scale.ones();
  }
}

inline
std::pair<double, arma::vec>
unstandardize(double&& intercept,
              arma::vec&& beta,
              const arma::rowvec& X_center,
              const arma::rowvec& X_scale,
              const double y_center,
              const double y_scale,
              const bool fit_intercept)
{
  using namespace arma;

  uword p = beta.n_rows;

  for (decltype(p) j = 0; j < p; ++j)
    beta(j) *= y_scale/X_scale(j);

  double X_bar_beta_sum = as_scalar(X_center * beta);

  if (fit_intercept)
    intercept = intercept*y_scale + y_center - X_bar_beta_sum;

  return std::make_pair(std::move(intercept),
                        std::move(beta));
}

#endif /* GOLEM_UTILS */
