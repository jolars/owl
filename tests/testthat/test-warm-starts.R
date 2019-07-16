test_that("warm starts reduce number of passes if used properly", {
  set.seed(3)

  xy <- golem:::randomProblem(1000, 10)
  x <- xy$x
  y <- xy$y

  model <- golem(penalty = "slope", sigma = 1)
  model$fit(x[1:500, ], y[1:500])
  passes_1 <- model$passes

  model$fit(x[501:1000, ], y[501:1000], warm_start = TRUE, sigma = 1)
  passes_2 <- model$passes

  expect_gt(passes_1, passes_2)
})
