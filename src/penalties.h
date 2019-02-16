#ifndef GOLEM_PENALTY_
#define GOLEM_PENALTY_

class Penalty {
public:
  virtual arma::vec eval(arma::vec y, const double L) = 0;

  virtual double loss(const arma::vec& beta) = 0;
};

class SLOPE : public Penalty {
public:
  SLOPE(const arma::vec& lambda) : lambda(lambda) {};

  arma::vec
  eval(arma::vec beta, const double L)
  {
    using namespace arma;

    uword p = beta.n_elem;

    // collect sign of beta and work with sorted absolutes
    vec beta_sign = sign(beta);
    beta = abs(beta);
    uvec beta_order = stable_sort_index(beta, "descend");
    beta = (beta(beta_order)).eval();

    vec s(p);
    vec w(p);
    vec betax(p);

    uvec idx_i(p);
    uvec idx_j(p);

    uword k = 0;

    for (uword i = 0; i < p; i++) {
      idx_i(k) = i;
      idx_j(k) = i;
      s(k)     = beta(i) - lambda(i)/L;
      w(k)     = s(k);

      while ((k > 0) && (w[k-1] <= w(k))) {
        k--;
        idx_j(k)  = i;
        s(k)     += s(k+1);
        w(k)      = s(k) / (i - idx_i(k) + 1.0);
      }
      k++;
    }

    for (uword j = 0; j < k; j++) {
      double d = std::max(w(j), 0.0);
      for (uword i = idx_i(j); i <= idx_j(j); i++) {
        betax(i) = d;
      }
    }

    // reset order
    betax(beta_order) = betax;

    // reset sign and return
    return (betax % beta_sign).eval();
  };

  double
  loss(const arma::vec& beta)
  {
    return arma::dot(lambda, arma::sort(arma::abs(beta), "descending"));
  }

private:
  const arma::vec& lambda;
};


// helper to choose prox
inline
std::unique_ptr<Penalty>
setupPenalty(const std::string& penalty_choice,
             const arma::vec& lambda)
{
  return std::unique_ptr<SLOPE>(new SLOPE{lambda});
}

#endif /* GOLEM_PENALTIES_ */
