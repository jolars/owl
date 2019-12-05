test_that("incorrect data dimensions throw errors", {
    expect_error(owl(matrix(1, 3, 3), double(2)))
})

test_that("na values in input throws errors", {
  expect_error(owl(matrix(NA, 3, 3), double(3)))
  expect_error(owl(matrix(3, 3, 3), c(NA, NA, 1)))
  expect_error(owl(matrix(3, 3, 3), double(3),
                   penalty = "group_slope",
                   groups = c(NA, 1, 1)))
})

test_that("orthogonalization and sparse input throws errors", {
  x <- Matrix::rsparsematrix(10, 10, 0.1)
  y <- double(10)
  groups <- sample(1:2, 10, replace = TRUE)
  expect_error(owl(x, y, groups = groups, penalty = "group_slope"))
})

test_that("erroneous group slope input throws errors", {
  x <- matrix(1, 3, 3)
  y <- double(3)
  groups <- c(1, 2, 2)

  expect_error(owl(x, y, groups = groups, penalty = "group_slope", sigma = -1))
  expect_error(owl(x, y, groups = groups, penalty = "group_slope",
                   lambda = 1:2, sigma = 1))
  expect_error(owl(x, y, groups = groups, penalty = "group_slope",
                   lambda = c(3, 1, -1)))
  expect_error(owl(x, y, groups = groups, penalty = "group_slope",
                   lambda = 1:3))
})

test_that("erroneous slope input throws errors", {
  x <- matrix(1, 3, 3)
  y <- double(3)

  expect_error(owl(x, y, penalty = "slope", lambda = 1:2))
  expect_error(owl(x, y, penalty = "slope", lambda = -c(1, 2)))
  expect_error(owl(x, y, penalty = "slope", lambda = 4:1))
})
