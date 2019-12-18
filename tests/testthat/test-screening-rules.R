test_that("screening rules return correct results for instances with known violations", {
  set.seed(216) # there is a violation for this seed and this setup


  for (family in c("gaussian", "binomial", "poisson", "multinomial")) {
    d <- owl:::randomProblem(100, 10, q = 0.1, response = family)

    coefs <- coef(owl::owl(d$x, d$y, family = family, screening = FALSE))

    y_mat <- matrix(NA, length(d$y), 3)
    y_mat[, 1] <- as.integer(d$y == 1)
    y_mat[, 2] <- as.integer(d$y == 2)
    y_mat[, 3] <- as.integer(d$y == 3)

    fit <- owl::owl(d$x, d$y, family = family, screening = TRUE)

    expect_equivalent(coefs, coef(fit), 1e-5)
  }
})

