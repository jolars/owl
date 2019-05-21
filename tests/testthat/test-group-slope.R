test_that("results from group slope mirror those from grpSLOPE package", {
  set.seed(1)
  x   <- matrix(rnorm(100^2), 100, 100)
  grp <- rep(rep(1:20), each = 5)
  b   <- c(runif(20), rep(0, 80))
  # (i.e., groups 1, 2, 3, 4, are truly significant)
  y   <- x %*% b + rnorm(10)
  fdr <- 0.1 # target false discovery rate
  sigma  <- 1

  golem_fit <- golem(x, y,
                     penalty = GroupSlope(groups = grp,
                                          fdr = fdr, sigma = sigma))

  gslope_fit <- grpSLOPE::grpSLOPE(x, y, group = grp, fdr = fdr, sigma = sigma)

  expect_equivalent(coef(golem_fit), coef(gslope_fit, scaled = FALSE),
                    tol = 1e-6)
})
