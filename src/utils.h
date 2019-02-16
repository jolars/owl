#ifndef GOLEM_UTILS
#define GOLEM_UTILS

template <typename T>
double
squaredNorm(const T& a)
{
  auto squared_norm = 0.0;

  for (const auto& a_i : a)
    squared_norm += a_i*a_i;

  return squared_norm;
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

#endif /* GOLEM_UTILS */
