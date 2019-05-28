test_that("results from group slope mirror those from grpSLOPE package", {
  set.seed(1)
  x   <- matrix(rnorm(100^2), 100, 100)
  grp <- rep(rep(1:20), each = 5)
  b   <- c(runif(20), rep(0, 80))
  # (i.e., groups 1, 2, 3, 4, are truly significant)
  y   <- x %*% b + rnorm(100)
  fdr <- 0.1 # target false discovery rate
  sigma  <- 1

  golem_fit <- golem(x, y,
                     penalty = "group_slope",
                     groups = grp,
                     fdr = fdr,
                     sigma = sigma)

  gslope_fit <- grpSLOPE::grpSLOPE(x, y, group = grp, fdr = fdr, sigma = sigma)

  expect_equivalent(coef(golem_fit), coef(gslope_fit, scaled = FALSE),
                    tol = 1e-5)
})

test_that("uneven group input is handled correctly", {
  set.seed(1)
  x   <- matrix(rnorm(100^2), 100, 100)
  grp <- sample(rep(rep(1:20), each = 5))
  b   <- ifelse(grp %in% 1:5, runif(sum(grp %in% 1:5)), 0)
  # (i.e., groups 1, 2, 3, 4, are truly significant)
  y   <- x %*% b + rnorm(nrow(x))
  fdr <- 0.1 # target false discovery rate
  sigma  <- 1

  golem_fit <- golem(x, y,
                     penalty = "group_slope",
                     groups = grp,
                     fdr = fdr,
                     sigma = sigma)

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
                     penalty = "group_slope",
                     groups = c(1, 1, 2),
                     sigma = 3)

  expect_equivalent(coef(golem_fit) != 0, c(TRUE, TRUE, TRUE, FALSE))
})

test_that("group_slope lambda sequences are computed properly", {
  set.seed(1)

  n <- 100
  p <- 5
  fdr <- 0.2

  x <- matrix(rnorm(p*n), n)
  y <- rnorm(n)

  set.seed(1)

  for (lambda in c("corrected", "mean", "max")) {
    groups <- sample(1:p, replace = TRUE)

    golem_fit <- golem(x, y, penalty = "group_slope", lambda = lambda,
                       fdr = fdr,
                       groups = groups)

    grps_lambda <- grpSLOPE:::lambdaGroupSLOPE(fdr = fdr, group = groups,
                                               wt = golem_fit@penalty@wt,
                                               n.obs = n,
                                               method = lambda)
    expect_equal(golem_fit@penalty@lambda,
                 grps_lambda)
  }



})

test_that("sigma estimation for Group SLOPE", {
  set.seed(1)

  problem <- golem:::random_problem(10, 5, sigma = 1)

  x <- problem$x
  y <- problem$y
  groups <- c(1, 1, 2, 2, 2)

  slope_fit <- grpSLOPE::grpSLOPE(x, y, group = groups, fdr = 0.2)
  golem_fit <- golem::golem(x, y, penalty = "group_slope",
                            fdr = 0.2,
                            groups = groups)

  expect_equal(slope_fit$sigma[length(slope_fit$sigma)],
               golem_fit@penalty@sigma)
})
