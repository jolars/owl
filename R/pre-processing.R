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
