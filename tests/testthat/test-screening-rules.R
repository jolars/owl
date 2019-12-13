test_that("screening rules return correct results for instances with known violations", {
  set.seed(216) # there is a violation for this seed and this setup
  d <- owl:::randomProblem(100, 10, q = 0.3)

  coefs <- coef(owl::owl(d$x, d$y, screening_rule = "none"))

  fit <- owl::owl(d$x, d$y, screening_rule = "strong")

  expect_equivalent(coefs, coef(fit), 1e-5)
})

