#ifndef GOLEM_PROX_
#define GOLEM_PROX_

using namespace arma;

class Prox {
public:
  virtual arma::vec eval(arma::vec y, double L) = 0;
};

class SLOPE : public Prox {
public:
  SLOPE(const arma::vec& lambda) : lambda(lambda) {};

  arma::vec
  eval(arma::vec y, double L)
  {
    using namespace arma;

    uword n = y.n_elem;

    // collect sign of y and work with sorted absolutes
    vec y_sign = sign(y);
    y = abs(y);
    uvec y_order = stable_sort_index(y, "descend");
    y = (y(y_order)).eval();

    vec s(n);
    vec w(n);
    vec x(n);

    uvec idx_i(n);
    uvec idx_j(n);

    sword k = 0;

    for (uword i = 0; i < n; i++) {
      idx_i(k) = i;
      idx_j(k) = i;
      s(k)     = y(i) - lambda(i)/L;
      w(k)     = s(k);

      while ((k > 0) && (w[k-1] <= w(k))) {
        k--;
        idx_j(k)  = i;
        s(k)     += s(k+1);
        w(k)      = s(k) / static_cast<double>(i - idx_i(k) + 1);
      }
      k++;
    }

    for (uword j = 0; j < k; j++) {
      double d = w(j);
      if (d < 0.0)
        d = 0.0;
      for (uword i = idx_i(j); i <= idx_j(j); i++) {
        x(i) = d;
      }
    }

    x(y_order) = x;

    // return sign and order
    return (x % y_sign).eval();
  };
private:
  arma::vec lambda;
};

#endif /* GOLEM_PROX_ */
