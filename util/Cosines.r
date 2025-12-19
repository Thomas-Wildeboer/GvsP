library(purrr)

# Cosine similarity between leading eigenvectors (PC1) of covariance matrices
Cosines <- function(a_list, g_list, p_list) {

  # leading eigenvector, normalized to unit length
  lead_vec <- function(M) {
    ev <- eigen(M, symmetric = TRUE)
    v <- ev$vectors[, 1]
    v / sqrt(sum(v^2))
  }

  # cosine similarity between two leading eigenvectors
  # abs() removes the arbitrary sign flip of eigenvectors
  cos_lead <- function(M1, M2) {
    v1 <- lead_vec(M1)
    v2 <- lead_vec(M2)
    abs(sum(v1 * v2))
  }

  p_a <- map2_dbl(p_list, a_list, cos_lead)
  g_a <- map2_dbl(g_list, a_list, cos_lead)
  p_g <- map2_dbl(p_list, g_list, cos_lead)

  out <- data.frame(
    Cosine = c(mean(p_a), mean(g_a), mean(p_g)),
    SE     = c(sd(p_a)/sqrt(length(p_a)),
               sd(g_a)/sqrt(length(g_a)),
               sd(p_g)/sqrt(length(p_g)))
  )
  rownames(out) <- c("P/A", "G/A", "P/G")
  
  return(out)
}
