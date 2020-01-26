test_that("SLOPE and owl agree for gaussian designs", {
  set.seed(1)

  problem <- owl:::randomProblem(1000, 50, sigma = 1)

  x <- scale(problem$x)
  y <- problem$y

  set.seed(0)
  res_slope <- SLOPE::SLOPE(x, y, sigma = 1, normalize = FALSE,
                            fdr = 0.25)
  set.seed(0)
  g <- owl(x, y,
           intercept = FALSE,
           sigma = 1,
           lambda = "gaussian",
           diagnostics = TRUE,
           standardize_features = FALSE,
           q = 0.25)

  expect_equivalent(res_slope$beta, coef(g), tol = 1e-5)
})

test_that("SLOPE and owl agree when computing lambda sequences", {
  set.seed(1)

  problem <- owl:::randomProblem(10, 5, sigma = 1)

  x <- problem$x
  y <- problem$y

  for (lambda in c("bh", "gaussian")) {
    slope_lambda <- SLOPE::SLOPE(x, y, sigma = 1, lambda = lambda,
                                 fdr = 0.1)$lambda
    owl_lambda <- owl(x, y, sigma = 1, lambda = lambda,
                      q = 0.1)$lambda*nrow(x)
    expect_equivalent(owl_lambda, slope_lambda)
  }
})
