preprocessResponse <- function(object, ...) {
  UseMethod("preprocessResponse", object)
}

preprocessResponse.Gaussian <- function(object, y) {
  y <- as.numeric(y)

  if (NCOL(y) > 1)
    stop("response for Gaussian regression must be one-dimensional.")

  y_center <- mean(y)
  y_scale  <- 1

  y <- as.matrix(y - y_center)

  list(y, y_center, y_scale, n_classes = 1L, class_names = NA_character_)
}

preprocessResponse.Binomial <- function(object, y) {
  if (NCOL(y) > 1)
    stop("response for binomial regression must be one-dimensional.")

  if (length(unique(y)) > 2)
    stop("more than two classes in response")

  if (length(unique(y)) == 1)
    stop("only one class in response.")

  y_table <- table(y)
  min_class <- min(y_table)

  if (min_class <= 1)
    stop("one class only has ", min_class, " observations.")

  class_names <- names(y_table)

  # Transform response to {-1, 1}, which is used internally
  y <- as.matrix(ifelse(as.numeric(as.factor(y)) == 1, -1, 1))

  list(y,
       y_center = 0,
       y_scale = 1,
       n_classes = 1L,
       class_names = class_names)
}
