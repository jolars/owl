Gaussian <- function() {
  structure(list(name = "gaussian"),
            class = c("Gaussian", "Family"))
}

Binomial <- function() {
  structure(list(name = "binomial"),
            class = c("Binomial", "Family"))
}

Poisson <- function() {
  structure(list(name = "poisson"),
            class = c("Poisson", "Family"))
}
