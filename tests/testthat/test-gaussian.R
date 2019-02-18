context("test-gaussian")

test_that("multiplication works", {
  set.seed(1)
  x1 <- rnorm(1000)
  x2 <- rnorm(1000)
  y <- 1 + 2*x1 + 3*x2

  df <- data.frame(y = y, x1 = x1, x2 = x2)
  lm_fit <- lm(y ~ x1 + x2)
  golem_fit <- golem::golem(cbind(x1, x2), y, family = "gaussian", lambda = rep(0, 2), intercept = TRUE)

  expect_equivalent(coef(lm_fit),
                    coef(golem_fit),
                    tol = 0.05)
})
