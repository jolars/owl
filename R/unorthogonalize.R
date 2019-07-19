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
