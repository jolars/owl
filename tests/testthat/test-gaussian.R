test_that("unregularized gaussian models work as expected", {
  set.seed(1)

  x <- as.matrix(abalone$x)
  y <- abalone$y

  lm_fit <- lm(y ~ x)
  golem_fit <- golem::golem(x, y, family = "gaussian",
                            penalty = Slope(sigma = 0.01),
                            solver = Fista(tol = 1e-6))

  expect_equivalent(coef(lm_fit),
                    coef(golem_fit),
                    tol = 0.1)
})
