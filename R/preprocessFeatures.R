preprocessFeatures <- function(x, standardize) {
  x <- as.matrix(x)

  p <- NCOL(x)
  n <- NROW(x)

  x_center <- double(p)
  x_scale  <- rep.int(1, p)

  if (standardize %in% c("both", "features")) {
    x_center <- colMeans(x)
    x_scale  <- colNorms(x)

    x_scale[x_scale == 0] <- 1

    x <- sweep(x, 2, x_center, check.margin = FALSE)
    x <- sweep(x, 2, x_scale, "/", check.margin = FALSE)
  }

  attr(x, "center") <- x_center
  attr(x, "scale") <- x_scale
  x
}
