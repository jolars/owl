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

test_that("binomial group slope models work", {
  set.seed(1)
  library(golem)

  X <- scale(matrix(rnorm(3000), ncol = 3))
  x1 <- X[, 1]
  x2 <- X[, 2]
  x3 <- X[, 3]
  z <- 1 + 2*x1 + 3*x2 + x3
  pr <- 1/(1 + exp(-z))
  y <- rbinom(1000, 1, pr)

  df <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)
  golem_fit <- golem(cbind(x1, x2, x3), y,
                     family = "binomial",
                     penalty = GroupSlope(groups = c(1, 1, 2), sigma = 3))

  expect_equivalent(coef(golem_fit) != 0, c(TRUE, TRUE, TRUE, FALSE))
})
