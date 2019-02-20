context("prox")

test_that("prox_sorted_L1 agrees with isotone package", {
  if (!requireNamespace("isotone", quietly = TRUE))
    skip("isotone is not available")
  if (!requireNamespace("SLOPE", quietly = TRUE))
    skip("SLOPE is not available")

  prox_slope_isotone <- function(y, lambda) {
    n <- length(y)

    # Solve the quadratic programming problem:
    #   min ||y - lambda - x||_2 s.t. x_1 >= x_2 >= ... >= x_n
    result <- isotone::activeSet(cbind(2:n, 1:(n-1)),
                                 y = y - lambda,
                                 weights = rep(1, n))

    # Enforce non-negativity constraint.
    pmax(result$x, 0)
  }


  n <- 20
  mu <- 1.5*(n:1)
  y <- sort(abs(rnorm(n, mean = mu)), decreasing = TRUE)
  lambda <- n:1

  args <- list(lambda = lambda, n = n, p = n, fdr = 0.2, sigma = 1,
               sigma_type = "user", lambda_type = "user")

  expect_equivalent(golem:::prox_slope_cpp(y, args),
                    prox_slope_isotone(y, lambda))
  expect_equivalent(golem:::prox_slope_cpp(y, args),
                    SLOPE:::prox_sorted_L1(y, lambda))
})

