test_that("diagnostics are working properly", {
  xy <- golem:::randomProblem(100, 2, q = 1)

  fit <- golem(xy$x, xy$y, diagnostics = TRUE, n_sigma = 1, sigma = 1)

  expect_is(fit$diagnostics, "data.frame")
  expect_silent(dont_plot(plotDiagnostics(fit)))
})
