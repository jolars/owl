test_that("incorrect data dimensions throw errors", {
  g <- golem()

  expect_error(g$fit(matrix(1, 3, 3), double(2)))
})

test_that("na values in input throws errors", {
  g <- golem()

  expect_error(g$fit(matrix(NA, 3, 3), double(3)))
  expect_error(g$fit(matrix(3, 3, 3), c(NA, NA, 1)))
  expect_error(g$fit(matrix(3, 3, 3), double(3),
                     penalty = "group_slope",
                     groups = c(NA, 1, 1)))
})

test_that("orthogonalization and sparse input throws errors", {
  x <- Matrix::rsparsematrix(10, 10, 0.1)
  y <- double(10)
  groups <- sample(1:2, 10, replace = TRUE)
  g <- golem(penalty = "group_slope")
  expect_error(g$fit(x, y, groups = groups))
})

test_that("erroneous group slope input throws errors", {
  x <- matrix(1, 3, 3)
  y <- double(3)
  groups <- c(1, 2, 2)
  g <- golem(penalty = "group_slope")

  expect_error(g$fit(x, y, groups = groups, sigma = -1))
  expect_error(g$fit(x, y, groups = groups, lambda = 1:2, sigma = 1))
  expect_error(g$fit(x, y, groups = groups, lambda = c(3, 1, -1)))
  expect_error(g$fit(x, y, groups = groups, lambda = 1:3))
})

test_that("erroneous slope input throws errors", {
  x <- matrix(1, 3, 3)
  y <- double(3)
  g <- golem(penalty = "slope")

  expect_error(g$fit(x, y, lambda = 1:2))
  expect_error(g$fit(x, y, lambda = -c(1, 2)))
  expect_error(g$fit(x, y, lambda = 4:1))
})
