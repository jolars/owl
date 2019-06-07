test_that("unregularized gaussian models work as expected", {
  set.seed(1)

  x <- as.matrix(abalone$x)
  y <- abalone$y
  lm_fit <- lm(y ~ x)
  g <- golem(family = "gaussian",
             sigma = 1e-4,
             penalty = "slope")
  g$fit(x, y)

  expect_equivalent(coef(lm_fit),
                    g$coef(),
                    tol = 1e-4)
})
