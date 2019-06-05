#' @include penalties.R

orthogonalizeGroups <- function(x,
                                group_id,
                                x_center,
                                x_scale,
                                standardize_features) {
  getGroupQR <- function(ids) {
    submat <- x[, ids, drop = FALSE]

    if (length(ids) == 1) {
      Q <- submat
      R <- 1
      P <- 1
    } else {

      if (inherits(submat, "sparseMatrix")) {
        submat_qr <- Matrix::qr(submat)
        Q <- Matrix::qr.Q(submat_qr)
        R <- Matrix::qrR(submat_qr)
        P <- submat_qr@q + 1
      } else {
        submat_qr <- qr(as.matrix(submat), LAPACK = TRUE)
        Q <- qr.Q(submat_qr)
        R <- qr.R(submat_qr)
        P <- submat_qr$pivot
      }
    }
    list(Q = Q, R = R, P = P)
  }

  lapply(group_id, getGroupQR)
}

orthogonalize <- function(x, penalty, x_center, x_scale, standardize_features) {

  orthogonalize <- penalty@orthogonalize
  group_id      <- penalty@group_id

  n_groups <- length(group_id)
  group_names <- names(group_id)

  n <- NROW(x)
  # p <- NCOL(x)

  if (orthogonalize) {
    ortho <- orthogonalizeGroups(x,
                                 group_id,
                                 x_center,
                                 x_scale,
                                 standardize_features)

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
    block_start <- utils::head(c(1, block_end + 1), n_groups)

    for (i in seq_len(n_groups)) {
      ind <- block_start[i]:block_end[i]
      grp[ind] <- group_names[i]
      x[, ind] <- as.matrix(ortho[[i]]$Q)
    }

    ortho_group_id <- groupID(grp)
    # set prior weights per group:
    wt <- sqrt(ortho_group_length)
    wt_per_coef <- rep(NA, ncol(x))

    for (i in 1:n_groups)
      wt_per_coef[ortho_group_id[[i]]] <- wt[i]

    penalty@ortho <- ortho
    penalty@ortho_group_length <- ortho_group_length
    penalty@ortho_group_id <- ortho_group_id

  } else {
    # set prior weights per group:
    wt <- sqrt(lengths(group_id))
    wt_per_coef <- rep(NA, ncol(x))

    for (i in seq_len(n_groups))
      wt_per_coef[group_id[[i]]] <- wt[i]
  }

  penalty@wt <- wt
  penalty@wt_per_coef <- wt_per_coef

  list(x = x, penalty = penalty)
}

standardize <- function(x, standardize_features) {
    p <- NCOL(x)

    if (standardize_features) {
      x_center <- Matrix::colMeans(x)

      if (!inherits(x, "sparseMatrix")) {
        x <- sweep(x, 2, x_center)
        x_scale <- apply(x, 2, norm, type = "2")
        x_scale[x_scale == 0] <- 1
        x <- sweep(x, 2, x_scale, "/")
      } else {
        x_scale <- standardizedSparseColNorms(x, x_center)
        x_scale[x_scale == 0] <- 1
        for (j in seq_len(p))
          x[, j] <- x[, j]/x_scale[j]
      }
    } else {
      x_center <- rep.int(0, p)
      x_scale  <- rep.int(1, p)
    }

    list(x, x_center, x_scale)
}


setGeneric(
  "preprocessFeatures",
  function(penalty, x, x_center, x_scale, standardize_features)
    standardGeneric("preprocessFeatures")
)

setMethod(
  "preprocessFeatures",
  "Penalty",
  function(penalty, x, x_center, x_scale, standardize_features) {
    list(penalty = penalty, x = x)
  }
)

setMethod(
  "preprocessFeatures",
  "GroupSlope",
  function(penalty, x, x_center, x_scale, standardize_features) {
    # orthogonalize (if required)
    res <- orthogonalize(x, penalty, x_center, x_scale, standardize_features)

    list(penalty = res$penalty, x = res$x)
  }
)
