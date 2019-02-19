context("test-slope")

test_that("SLOPE and golem agree for gaussian designs", {
  set.seed(1)

  problem <- golem:::random_problem(1000, 50, sigma = 1)

  X <- problem$X
  y <- problem$y

  set.seed(0)
  res_slope <- SLOPE::SLOPE(X, y, sigma = 1, normalize = FALSE)
  set.seed(0)
  res_golem <- golem::golem(X, y, sigma = 1, intercept = FALSE,
                            standardize = "none")

  expect_equivalent(res_slope$beta, res_golem$coefficients, tol = 0.05)
})

test_that("SLOPE and golem agree when computing lambda sequences", {
  set.seed(1)

  problem <- golem:::random_problem(10, 5, sigma = 1)

  X <- problem$X
  y <- problem$y

  for (lambda in c("bhq", "gaussian")) {
    slope_lambda <- SLOPE::SLOPE(X, y, sigma = 1, lambda = lambda)$lambda
    golem_lambda <- golem::golem(X, y, sigma = 1, lambda = lambda)$lambda
    expect_equivalent(golem_lambda, slope_lambda)
  }
})

# test_that("fdr is kept for binomial family and orthogonal design", {
#   set.seed(1)
#
#   n <- 1000
#   p <- n/2
#   q <- 0.1
#   m <- floor(q*p)
#
#   x <- matrix(rnorm(n*p), n, p)
#   i <- sample(p, m)
#   beta <- rep(0, p)
#   beta[i] <- sqrt(2*log(p))
#
#   z <- x %*% beta
#   pr <- 1/(1 + exp(-z))
#   y <- rbinom(n, 1, pr)
#
#   golem_unreg <- golem::golem(x, y, family = "binomial", lambda = "bhq", sigma = 0, intercept = FALSE)
#   sigma <- sd(as.numeric(predict(golem_unreg, x, "response")) - sigmoid(y))
#
#   golem_fit <- golem::golem(x, y, family = "binomial", lambda = "bhq", sigma = sigma, intercept = FALSE)
#   selected <- golem_fit$coefficients > 0
#
#   lmfit <- lm(z ~ x)
#   sigma <- sqrt(sum(lmfit$residuals^2) / (n-p))
#
#   problem <- golem:::random_problem(5000, 100, 20, sigma = 1)
#
#   set.seed(1)
#
#   n <- 5000
#   p <- 100
#   q <- 0.2
#   m <- floor(q*p)
#   amplitude <- 3
#
#   x <- matrix(rnorm(n*p), n, p)
#   i <- sample(p, m)
#   beta <- amplitude * (1:p %in% i)
#   y <- x %*% beta + rnorm(n, sd = 1)
#
#   golem_fit <- golem::golem(X, y, family = "gaussian", lambda = "bhq", sigma = 1, intercept = FALSE)
#
#   slope_fit <- SLOPE::SLOPE(X, y, lambda = "bhq", fdr = 0.2)
#
#   selected <- slope_fit$selected
#   #selected <- which(coef(golem_fit) > 0)
#   V <- length(setdiff(selected, i))
#   R <- length(selected)
#   fdr <- V/R
#   fdr
#
#   Sx2 <- x[, selected]
#   f <- glm(y ~ x2, family = "binomial")
#
#   c(golem_fit$intercept, golem_fit$coefficients[1])
#   f$coefficients[1:2]
#
# })
