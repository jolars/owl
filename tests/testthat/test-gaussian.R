test_that("unregularized gaussian models work as expected", {
  set.seed(1)

  x <- as.matrix(abalone$x)
  y <- abalone$y
  lm_fit <- lm(y ~ x)
  g <- owl(x,
           y,
           family = "gaussian",
           sigma = 1e-9)

  expect_equivalent(coef(lm_fit),
                    coef(g),
                    tol = 1e-4)
})
