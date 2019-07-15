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

groupify <- function(x,
                     x_center,
                     x_scale,
                     groups,
                     standardize_features,
                     orthogonalize) {

  group_id <- groupID(groups)

  n_groups <- length(group_id)
  group_names <- names(group_id)

  n <- NROW(x)

  groups <- list()

  if (orthogonalize) {
    ortho <- orthogonalizeGroups(x,
                                 group_id,
                                 x_center,
                                 x_scale,
                                 standardize_features)

    # determine sizes of orthogonalized groups:
    ortho_group_length <- rep(NA, n_groups)
    names(ortho_group_length) <- group_names

    for (i in 1:n_groups)
      ortho_group_length[i] <- ncol(ortho[[i]]$Q)

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

    groups$ortho <- ortho
    groups$ortho_group_length <- ortho_group_length
    groups$ortho_group_id <- ortho_group_id

  } else {
    # set prior weights per group:
    wt <- sqrt(lengths(group_id))
    wt_per_coef <- rep(NA, ncol(x))

    for (i in seq_len(n_groups))
      wt_per_coef[group_id[[i]]] <- wt[i]
  }

  groups$wt <- wt
  groups$wt_per_coef <- wt_per_coef
  groups$group_id <- group_id
  groups$orthogonalize <- orthogonalize

  list(x, groups)
}


unorthogonalize <- function(betas,
                            groups) {

  group_id <- groups$group_id
  n_groups <- length(group_id)
  ortho_group_id <- groups$ortho_group_id

  group_lengths <- lengths(group_id)
  ortho_group_lengths <- lengths(ortho_group_id)

  ortho <- groups$ortho

  if (all(group_lengths == ortho_group_lengths)) {

    beta_tilde <- betas

    for (i in seq_len(n_groups)) {
      # c corresponds to the (reordered) group
      # structure in orthogonalized version of X
      ci <- betas[ortho_group_id[[i]]]
      li <- ortho_group_lengths[i]

      if (inherits(ortho[[i]]$R, "sparseMatrix")) {
        bi <- tryCatch({
          Matrix::solve(ortho[[i]]$R, ci)
        }, error = function(err) {
          warning(paste("golem caught an error:", err))
          rep(NA, li)
        })
      } else {
        bi <- tryCatch({
          backsolve(ortho[[i]]$R, ci)
        }, error = function(err) {
          warning(paste("golem caught an error:", err))
          rep(NA, li)
        })
      }

      or <- double(li)

      for (j in seq_len(li))
        or[j] <- which(ortho[[i]]$P == j)

      # beta corresponds to the group structure in the original matrix
      beta_tilde[group_id[[i]]] <- bi[or]
    }
  } else {
    beta_tilde <- NULL
  }

  beta_tilde
}
