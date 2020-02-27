test_that("screening rules return correct results for instances with known violations", {
  set.seed(216) # there is a violation for this seed and this setup

  for (family in c("gaussian", "binomial", "poisson", "multinomial")) {
    d <- owl:::randomProblem(100, 10, q = 0.1, response = family)

    fit0 <- owl(d$x, d$y, family = family, screening = FALSE)
    fit1 <- owl(d$x, d$y, family = family, screening = TRUE)

    expect_equivalent(coef(fit0), coef(fit1), 1e-4)
  }
})

test_that("basic screening rule works", {
  set.seed(213)

  xy <- owl:::randomProblem(100, 10, q = 0.1)

  fit <- owl(xy$x, xy$y, screening = TRUE)

  expect_lt(length(fit$active_sets[[1]]), ncol(xy$x))
})
