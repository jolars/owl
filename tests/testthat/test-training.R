test_that("model training works with trainGolem", {
  set.seed(48)

  for (family in c("gaussian", "binomial")) {
    xy <- golem:::randomProblem(1e2, 4, n_groups = 2, response = family)

    for (penalty in c("slope", "group_slope", "lasso")) {

      fit <- trainGolem(xy$x,
                        xy$y,
                        penalty = penalty,
                        family = family,
                        groups = xy$groups,
                        number = 2,
                        fdr = c(0.1, 0.2),
                        repeats = 2,
                        n_sigma = 3,
                        n_lambda = 3)
      expect_s3_class(fit, "TrainedGolem")
    }
  }
})

test_that("erroneous input throws errors in plot.trainGolem", {
  xy <- golem:::randomProblem(1e3, 2)
  x <- xy$x
  y <- xy$y

  fit <- trainGolem(x, y)

  expect_error(plot(fit, measure = "auc"))
  p <- plot(fit)
  expect_s3_class(p, "trellis")
  expect_silent(dont_plot(p))

  fit <- trainGolem(x, y, fdr = c(0.1, 0.2))

  p <- plot(fit)
  expect_s3_class(p, "trellis")
})
