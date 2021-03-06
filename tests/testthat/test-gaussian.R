test_that("unregularized gaussian models work as expected", {
  set.seed(1)

  x <- as.matrix(abalone$x)
  y <- abalone$y

  lm_fit <- lm(y ~ as.matrix(x))

  g <- owl(x,
           y,
           family = "gaussian",
           sigma = 1e-12)

  expect_equivalent(coef(lm_fit),
                    coef(g),
                    tol = 1e-3)
})

test_that("wide and tall inputs work correctly", {

  set.seed(926)

  grid <- expand.grid(n = c(50, 100),
                      p = c(50, 100),
                      density = c(1, 0.5))

  for (i in seq_len(nrow(grid))) {
    n <- grid$n[i]
    p <- grid$p[i]
    d <- grid$density[i]

    xy <- owl:::randomProblem(n, p, density = d)

    expect_silent(owl(xy$x, xy$y, n_sigma = 5))
  }
})
