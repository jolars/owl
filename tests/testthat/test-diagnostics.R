test_that("diagnostics are working properly", {
  xy <- golem:::randomProblem(100, 2, q = 1)

  model <- golem(diagnostics = TRUE, n_sigma = 1)
  model$fit(xy$x, xy$y)

  expect_is(model$diagnostics, "Diagnostics")
  expect_is(model$diagnostics$data, "data.frame")
  expect_silent(dont_plot(model$diagnostics$plot()))
})
