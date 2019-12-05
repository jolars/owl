test_that("screening rules return correct results for instances with known violations", {
  set.seed(216) # there is a violation for this seed and this setup
  d <- owl:::randomProblem(100, 10, q = 0.3)

  fit <- owl::owl(d$x, d$y, screening_rule = "none")
  coefs <- coef(fit)

  for (rule in c("none", "safe", "strong")) {
    fit <- owl::owl(d$x, d$y, screening_rule = rule, diagnostics = TRUE)
    expect_equivalent(coefs, coef(fit), 1e-5)
  }
})

