context("test-slope")

test_that("multiplication works", {
  if (!requireNamespace("SLOPE", quietly = TRUE))
    skip("SLOPE package not available")

  set.seed(1)

  problem <- golem:::random_problem(10000, 50, sigma = 1)

  X <- problem$X
  y <- problem$y

  set.seed(0)
  res_slope <- SLOPE::SLOPE(X, y, sigma = 1)
  set.seed(0)
  res_golem <- golem:::golem(X, y, sigma = 1, intercept = FALSE)

  expect_equivalent(res_slope$beta, res_golem$beta, tol = 0.05)
})

test_that("fdr is kept for binomial family and orthogonal design", {
  set.seed(1)

  n <- 10000
  p <- 100
  m <- 5

  x <- matrix(rnorm(n*p), n, p)
  i <- sample(p, m)

  z <- 1 + x[, 1]*2
  pr <- 1/(1 + exp(-z))
  y <- rbinom(n, 1, pr)

  golem_fit <- golem::golem(x, y, family = "binomial", lambda = "bhq", sigma = NULL)
  f <- glm(y ~ x, family = "binomial")

  c(golem_fit$intercept, golem_fit$coefficients[1])
  f$coefficients[1:2]

})


res <- svd(crossprod(X))
