test_that("different solvers return equivalent results", {
  xy <- owl:::randomProblem()

  solvers <- c("admm", "fista")

  set.seed(-1)

  for (i in seq_along(solvers)) {
    fit <- owl(xy$x, xy$y, intercept = FALSE, solver = solvers[i], sigma = 1)

    if (i != 1)
      expect_equivalent(coef(fit), coef(fit_old), tol = 1e-5)

    fit_old <- fit
  }
})
