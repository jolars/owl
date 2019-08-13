test_that("sparse and dense implementations give equivalent results", {
  set.seed(2)
  n <- 1000
  p <- 2

  for (family in c("gaussian", "binomial")) {
    for (standardize in c(TRUE, FALSE)) {
      d <- prague:::randomProblem(n, p, 0.5, density = 0.5, response = family)
      x <- d$x
      y <- d$y
      beta <- d$beta

      g <- golem(x, y, family = family, standardize_features = standardize)

      sparse_coefs <- coef(g)

      g <- golem(as.matrix(x), y, family = family,
                 standardize_features = standardize)

      dense_coefs <- coef(g)

      expect_equal(sparse_coefs, dense_coefs, tol = 1e-4)
    }
  }

  p <- 10
  d <- prague:::randomProblem(n, p, 0.5, density = 0.5, response = "gaussian",
                             n_groups = 5)
  x <- d$x
  y <- d$y
  beta <- d$beta

  for (orthogonalize in c(TRUE, FALSE)) {
    sparse_coefs <- coef(golem(x, y,
                               sigma = 1,
                               penalty = "group_slope",
                               groups = d$groups,
                               standardize_features = FALSE,
                               orthogonalize = orthogonalize))
    dense_coefs <- coef(golem(as.matrix(x), y,
                              sigma = 1,
                              penalty = "group_slope",
                              groups = d$groups,
                              standardize_features = FALSE,
                              orthogonalize = orthogonalize))
    expect_equal(sparse_coefs, dense_coefs, tol = 1e-7)
  }
})
