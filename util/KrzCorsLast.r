# Krzanowski common subspaces for last n/2-1 dimensions

KrzCorsLast <- function (a_list, g_list, p_list) {
  full.dim <- dim(a_list[[1]])[1]
  ret.dim <- round((full.dim/2)-2)

  cors_p_a <- list()
  cors_g_a <- list()
  cors_p_g <- list()

  for (i in 1:length(a_list)) {
    cov.a <- a_list[[i]]
    cov.p <- p_list[[i]]
    cov.g <- g_list[[i]]

    eVec.a <- eigen(cov.a)$vectors
    eVec.p <- eigen(cov.p)$vectors
    eVec.g <- eigen(cov.g)$vectors

    # Krzanowski common subspaces (S = A^T B B^T A)
    S.p <- t(eVec.a[,(full.dim-ret.dim):full.dim]) %*%
            eVec.p[,(full.dim-ret.dim):full.dim] %*%
            t(eVec.p[,(full.dim-ret.dim):full.dim]) %*%
            eVec.a[,(full.dim-ret.dim):full.dim]

    S.g <- t(eVec.a[,(full.dim-ret.dim):full.dim]) %*%
            eVec.g[,(full.dim-ret.dim):full.dim] %*%
            t(eVec.g[,(full.dim-ret.dim):full.dim]) %*%
            eVec.a[,(full.dim-ret.dim):full.dim]

    S.pg <- t(eVec.p[,(full.dim-ret.dim):full.dim]) %*%
            eVec.g[,(full.dim-ret.dim):full.dim] %*%
            t(eVec.g[,(full.dim-ret.dim):full.dim]) %*%
            eVec.p[,(full.dim-ret.dim):full.dim]

    # sum eigenvalues of S and divide by number of dimensions
    cors_p_a[[i]] <- sum(eigen(S.p)$values)/(ret.dim+1)
    cors_g_a[[i]] <- sum(eigen(S.g)$values)/(ret.dim+1)
    cors_p_g[[i]] <- sum(eigen(S.pg)$values)/(ret.dim+1)
  }

  cors_g_se <- sd(unlist(cors_g_a))/sqrt(length(unlist(cors_g_a)))
  cors_p_se <- sd(unlist(cors_p_a))/sqrt(length(unlist(cors_p_a)))
  cors_pg_se <- sd(unlist(cors_p_g))/sqrt(length(unlist(cors_p_g)))

  krzcors <- data.frame(
      c(mean(unlist(cors_p_a)), mean(unlist(cors_g_a)), mean(unlist(cors_p_g))),
      c(cors_p_se, cors_g_se, cors_pg_se)
  )

  # set column names
  colnames(krzcors) <- c(
      "KrzCor Last",
      "SE"
  )

  # set row names
  rownames(krzcors) <- c(
      "P,A",
      "G,A",
      "P,G"
  )

  return(krzcors)
}