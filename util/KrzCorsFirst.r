# Krzanowski common subspaces for first n/2-1 dimensions

KrzCorsFirst <- function (a_list, g_list, p_list) {
  full.dim <- dim(a_list[[1]])[1]
  ret.dim <- round((full.dim/2)-1)

  cors_p_a <- list()
  cors_g_a <- list()
  cors_g_p <- list()

  for (i in 1:length(a_list)) {
    cov.a <- a_list[[i]]
    cov.p <- p_list[[i]]
    cov.g <- g_list[[i]]

    eVec.a <- eigen(cov.a)$vectors
    eVec.p <- eigen(cov.p)$vectors
    eVec.g <- eigen(cov.g)$vectors

    # Krzanowski common subspaces (S = A^T B B^T A)
    S.p <- t(eVec.a[,1:ret.dim]) %*%
            eVec.p[,1:ret.dim] %*%
            t(eVec.p[,1:ret.dim]) %*%
            eVec.a[,1:ret.dim]

    S.g <- t(eVec.a[,1:ret.dim]) %*%
            eVec.g[,1:ret.dim] %*%
            t(eVec.g[,1:ret.dim]) %*%
            eVec.a[,1:ret.dim]

    S.pg <- t(eVec.p[,1:ret.dim]) %*%
            eVec.g[,1:ret.dim] %*%
            t(eVec.g[,1:ret.dim]) %*%
            eVec.p[,1:ret.dim]

    # sum eigenvalues of S and divide by number of dimensions
    cors_p_a[[i]] <- sum(eigen(S.p)$values)/ret.dim
    cors_g_a[[i]] <- sum(eigen(S.g)$values)/ret.dim
    cors_g_p[[i]] <- sum(eigen(S.pg)$values)/ret.dim
  }

  cors_g_se <- sd(unlist(cors_g_a))/sqrt(length(unlist(cors_g_a)))
  cors_p_se <- sd(unlist(cors_p_a))/sqrt(length(unlist(cors_p_a)))
  cors_gp_se <- sd(unlist(cors_g_p))/sqrt(length(unlist(cors_g_p)))

  krzcors <- data.frame(
      c(mean(unlist(cors_p_a)), mean(unlist(cors_g_a)), mean(unlist(cors_g_p))),
      c(cors_p_se, cors_g_se, cors_gp_se)
  )

  # set column names
  colnames(krzcors) <- c(
      "KrzCor First",
      "SE"
  )

  # set row names
  rownames(krzcors) <- c(
      "P,A",
      "G,A",
      "G,P"
  )

  return(krzcors)
}