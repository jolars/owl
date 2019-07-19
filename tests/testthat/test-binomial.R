test_that("unregularized logistic regression matches output from glm()", {
  set.seed(1)
  X <- scale(matrix(rnorm(3000), ncol = 3))
  x1 <- X[, 1]
  x2 <- X[, 2]
  x3 <- X[, 3]
  z <- 1 + 2*x1 + 3*x2 + x3
  pr <- 1/(1 + exp(-z))
  y <- rbinom(1000, 1, pr)

  df <- data.frame(y = y, x1 = x1, x2 = x2)
  glm_fit <- glm(y ~ x1 + x2 + x3, data = df, family = "binomial")

  g_model <- golem(cbind(x1, x2, x3), y,
                   family = "binomial",
                   penalty = "slope",
                   diagnostics = TRUE,
                   sigma = 0)

  expect_equivalent(coef(glm_fit),
                    coef(g_model),
                    tol = 0.01)
})

test_that("unregularized group slope logistic regression matches output from glm()", {
  set.seed(1)
  X <- scale(matrix(rnorm(3000), ncol = 3))
  x1 <- X[, 1]
  x2 <- X[, 2]
  x3 <- X[, 3]
  z <- 1 + 2*x1 + 3*x2 + x3
  pr <- 1/(1 + exp(-z))
  y <- rbinom(1000, 1, pr)

  df <- data.frame(y = y, x1 = x1, x2 = x2)
  glm_fit <- glm(y ~ x1 + x2 + x3, data = df, family = "binomial")
  gol_fit <- golem(cbind(x1, x2, x3), y,
                   groups = c(1, 1, 2),
                   family = "binomial",
                   penalty = "group_slope",
                   diagnostics = TRUE,
                   sigma = 0)

  expect_equivalent(coef(glm_fit),
                    coef(gol_fit),
                    tol = 0.01)
})
