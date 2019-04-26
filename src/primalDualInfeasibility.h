#pragma once

#include <RcppArmadillo.h>
#include "families.h"
#include "penalties.h"
#
# inline
# void
# primalDualInfeasibility(const std::unique_ptr<Gaussian>& family,
#                         const std::unique_ptr<Lasso>& penalty,
#                         const arma::mat& lin_pred,
#                         const arma::mat& y,
#                         const arma::mat& beta,
#                         const arma::mat& grad,
#                         double& primal,
#                         double& dual,
#                         double& infeasibility)
# {
#   double f = family->loss(lin_pred, y);
#
#   primal = f + penalty->lambda*arma::norm(beta, 1);
#   dual =
# }
