#' @include penalties.R

orthogonalizeGroups <- function(x, group_id) {
  getGroupQR <- function(ids) {
    submat <- as.matrix(x[, ids, drop = FALSE])

    if (length(ids) == 1) {
      Q <- submat
      R <- 1
      P <- 1
    } else {
      submat_qr <- qr(submat, LAPACK = TRUE)
      Q <- qr.Q(submat_qr)
      R <- qr.R(submat_qr)
      P <- submat_qr$pivot
    }
    list(Q = Q, R = R, P = P)
  }

  lapply(group_id, getGroupQR)
}

orthogonalize <- function(x, penalty) {

  orthogonalize <- penalty@orthogonalize
  group_id      <- penalty@group_id

  n_groups <- length(group_id)
  group_names <- names(group_id)

  n <- NROW(x)
  p <- NCOL(x)
  x_attr <- attributes(x)

  if (orthogonalize) {
    ortho <- orthogonalizeGroups(x, group_id)

    # determine sizes of orthogonalized groups:
    ortho_group_length <- rep(NA, n_groups)
    names(ortho_group_length) <- group_names

    for (i in 1:n_groups) {
      ortho_group_length[i] <- ncol(ortho[[i]]$Q)
    }

    # overwrite x with a matrix that contains orthogonalizations
    # of the original blocks of x
    # and set up grouping info for the orthogonalized version of x
    # to be used in the optimization
    x <- matrix(nrow = n, ncol = sum(ortho_group_length))
    grp <- rep(NA, sum(ortho_group_length))
    block_end <- cumsum(ortho_group_length)
    block_start <- head(c(1, block_end + 1), n_groups)

    for (i in seq_len(n_groups)) {
      ind <- block_start[i]:block_end[i]
      grp[ind] <- group_names[i]
      x[, ind] <- ortho[[i]]$Q
    }

    ortho_group_id <- groupID(grp)
    # set prior weights per group:
    wt <- sqrt(ortho_group_length)
    wt_per_coef <- rep(NA, ncol(x))

    for (i in 1:n_groups) {
      wt_per_coef[ortho_group_id[[i]]] <- wt[i]
    }

    penalty@ortho <- ortho
    penalty@ortho_group_length <- ortho_group_length

    attr(x, "center") <- x_attr$center
    attr(x, "scale") <- x_attr$scale

  } else {
    # set prior weights per group:
    wt <- sqrt(lengths(group_id))
    wt_per_coef <- rep(NA, ncol(x))

    for (i in seq_len(n_groups)) {
      wt_per_coef[group_id[[i]]] <- wt[i]
    }
  }

  penalty@wt <- wt
  penalty@ortho_group_id <- ortho_group_id
  penalty@wt_per_coef <- wt_per_coef

  attr(x, "n") <- n
  attr(x, "p") <- p

  list(x = x, penalty = penalty)
}

standardize <- function(x, standardize) {
    x <- as.matrix(x)

    p <- NCOL(x)
    n <- NROW(x)

    x_center <- double(p)
    x_scale  <- rep.int(1, p)

    if (standardize %in% c("both", "features")) {
      x_center <- colMeans(x)

      x <- sweep(x, 2, x_center, check.margin = FALSE)

      x_scale  <- colNorms(x)
      x_scale[x_scale == 0] <- 1

      x <- sweep(x, 2, x_scale, "/", check.margin = FALSE)
    }

    attr(x, "n")      <- n
    attr(x, "p")      <- p
    attr(x, "center") <- x_center
    attr(x, "scale")  <- x_scale

    x
}


setGeneric("preprocessFeatures",
           function(penalty, x, standardize)
             standardGeneric("preprocessFeatures"))

setMethod(
  "preprocessFeatures",
  "Penalty",
  function(penalty, x, standardize) {
    x <- standardize(x, standardize)
    list(x = x, penalty = penalty)
  }
)

setMethod(
  "preprocessFeatures",
  "GroupSlope",
  function(penalty, x, standardize) {
    # first perform standardization if required
    x <- standardize(x, standardize)

    # orthogonalize (if required)
    res <- orthogonalize(x, penalty)

    list(x = res$x, penalty = res$penalty)
  }
)
