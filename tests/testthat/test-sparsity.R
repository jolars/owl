test_that("sparse and dense implementations give equivalent results", {
  set.seed(2)
  n <- 1000
  p <- 2

  g <- golem(sigma = 1)

  for (family in c("gaussian", "binomial")) {
    for (standardize in c(TRUE, FALSE)) {
      d <- golem:::randomProblem(n, p, 0.5, density = 0.5, response = family)
      x <- d$x
      y <- d$y
      beta <- d$beta

      g$fit(x, y, family = family, standardize_features = standardize)
      sparse_coefs <- g$coef()
      g$fit(as.matrix(x), y)
      dense_coefs <- g$coef()

      expect_equal(sparse_coefs, dense_coefs, tol = 1e-4)
    }
  }

  p <- 10
  d <- golem:::randomProblem(n, p, 0.5, density = 0.5, response = "gaussian",
                             n_groups = 5)
  x <- d$x
  y <- d$y
  beta <- d$beta

  g <- golem(sigma = 1, penalty = "group_slope",
             standardize_features = FALSE)

  for (orthogonalize in c(TRUE, FALSE)) {
    sparse_coefs <- g$fit(x, y,
                          groups = d$groups,
                          orthogonalize = orthogonalize)$coef()
    dense_coefs <- g$fit(as.matrix(x), y,
                         groups = d$groups,
                         orthogonalize = orthogonalize)$coef()
    expect_equal(sparse_coefs, dense_coefs, tol = 1e-7)
  }
})
