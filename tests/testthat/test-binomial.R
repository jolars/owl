test_that("unregularized logistic regression matches output from glm()", {
  set.seed(1)
  X <- scale(matrix(rnorm(2000), ncol = 2))
  x1 <- X[, 1]
  x2 <- X[, 2]
  z <- 1 + 2*x1 + 3*x2
  pr <- 1/(1 + exp(-z))
  y <- rbinom(1000, 1, pr)

  df <- data.frame(y = y, x1 = x1, x2 = x2)
  glm_fit <- glm(y ~ x1 + x2, data = df, family = "binomial")
  golem_fit <- golem::golem(cbind(x1, x2), y, family = "binomial",
                            penalty = Slope(sigma = 0))

  expect_equivalent(coef(glm_fit),
                    coef(golem_fit),
                    tol = 0.01)
})
