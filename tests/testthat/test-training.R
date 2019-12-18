test_that("model training works with trainOwl", {
  set.seed(48)

  for (family in c("gaussian", "binomial")) {
    xy <- owl:::randomProblem(1e3, 4, n_groups = 2, response = family)

    for (penalty in c("slope", "group_slope")) {

      fit <- trainOwl(xy$x,
                      xy$y,
                      penalty = penalty,
                      family = family,
                      groups = xy$groups,
                      number = 2,
                      fdr = c(0.1, 0.2),
                      repeats = 2,
                      n_sigma = 2)
      expect_s3_class(fit, "TrainedOwl")
    }
  }
})

test_that("erroneous input throws errors in plot.trainOwl", {
  xy <- owl:::randomProblem(1e3, 2)
  x <- xy$x
  y <- xy$y

  fit <- trainOwl(x, y)

  expect_error(plot(fit, measure = "auc"))
  p <- plot(fit)
  expect_s3_class(p, "trellis")
  expect_silent(dont_plot(p))

  fit <- trainOwl(x, y, fdr = c(0.1, 0.2), number = 2)

  p <- plot(fit)
  expect_s3_class(p, "trellis")
})
