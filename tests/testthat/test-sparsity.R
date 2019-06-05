test_that("sparse and dense implementations give equivalent results", {
  set.seed(2)
  n <- 1000
  p <- 2

  for (family in c("gaussian", "binomial")) {
    for (standardize in c(TRUE, FALSE)) {
      d <- golem:::randomProblem(n, p, 0.5, density = 0.5, response = family)
      x <- d$x
      y <- d$y
      beta <- d$beta

      sparse_fit <- golem(x, y, sigma = 1, standardize_features = standardize)
      dense_fit  <- golem(as.matrix(x), y, sigma = 1,
                          standardize_features = standardize)

      expect_equal(coef(sparse_fit),
                   coef(dense_fit),
                   tol = 1e-4)
    }
  }

  p <- 10
  d <- golem:::randomProblem(n, p, 0.5, density = 0.5, response = "gaussian",
                             n_groups = 5)
  x <- d$x
  y <- d$y
  beta <- d$beta

  for (orthogonalize in c(TRUE, FALSE)) {
    sparse_fit <- golem(x, y, sigma = 1, penalty = "group_slope",
                        standardize_features = FALSE,
                        groups = d$groups, orthogonalize = orthogonalize)

    dense_fit <- golem(as.matrix(x), y, sigma = 1, penalty = "group_slope",
                       standardize_features = FALSE,
                       groups = d$groups, orthogonalize = orthogonalize)
    expect_equal(coef(sparse_fit),
                 coef(dense_fit),
                 tol = 1e-7)
  }
})
