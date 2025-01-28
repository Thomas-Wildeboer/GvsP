# Inputs CSV file of matrices from matrixsim.r and optional heritability for scaling
# Outputs R list of A, G, and P matrices for convenience

`read_matrices` <- function(matrices_csv, scale_p=FALSE, h, cor=FALSE) {
  a_list <- list()
  g_list <- list()
  p_list <- list()

  for(i in 1:nrow(matrices_csv)){
    a <- strsplit(matrices_csv[i, "A"], split=" ")[[1]]
    g <- strsplit(matrices_csv[i, "G"], split=" ")[[1]]
    p <- strsplit(matrices_csv[i, "P"], split=" ")[[1]]

    dim <- sqrt(length(a))

    a <- matrix(as.numeric(a), nrow=dim, ncol=dim)
    g <- matrix(as.numeric(g), nrow=dim, ncol=dim)
    p <- matrix(as.numeric(p), nrow=dim, ncol=dim)
    
    if (cor == TRUE) {
      # convert to correlation matrices
      a <- cov2cor(a)
      g <- cov2cor(g)
      p <- cov2cor(p)
    }

    if (scale_p == TRUE) {
      # adjust p by h (Cheverud 1988)
      p <- h * p
    }

    a_list[[i]] <- a
    g_list[[i]] <- g
    p_list[[i]] <- p
  }

  return(list(a_list, g_list, p_list))
}
