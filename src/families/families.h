#pragma once

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

#include "family.h"
#include "gaussian.h"
#include "binomial.h"
#include "poisson.h"
#include "multinomial.h"

// helper to choose family
inline std::unique_ptr<Family> setupFamily(const std::string& family_choice)
{
  if (family_choice == "binomial")
    return std::unique_ptr<Binomial>(new Binomial);
  else if (family_choice == "poisson")
    return std::unique_ptr<Poisson>(new Poisson);
  else if (family_choice == "multinomial")
    return std::unique_ptr<Multinomial>(new Multinomial);
  else
    return std::unique_ptr<Gaussian>(new Gaussian);
}

