library(purrr)

# Acute angle (in degrees) between leading eigenvectors (PC1) of covariance matrices:
Angles <- function(a_list, g_list, p_list) {

  rad2deg <- function(x) x * 180 / pi

  # leading eigenvector, normalized to unit length
  lead_vec <- function(M) {
    ev <- eigen(M, symmetric = TRUE)
    v <- ev$vectors[, 1]
    v / sqrt(sum(v^2))
  }

  # acute angle between two leading eigenvectors (degrees)
  angle_lead_deg <- function(M1, M2) {
    v1 <- lead_vec(M1)
    v2 <- lead_vec(M2)

    dot <- abs(drop(crossprod(v1, v2)))
    dot <- max(-1, min(1, dot))

    rad2deg(acos(dot))  # in [0, 90]
  }

  p_a <- map2_dbl(p_list, a_list, angle_lead_deg)
  g_a <- map2_dbl(g_list, a_list, angle_lead_deg)
  p_g <- map2_dbl(p_list, g_list, angle_lead_deg)

  out <- data.frame(
    Angle_deg = c(mean(p_a), mean(g_a), mean(p_g)),
    SE_deg    = c(sd(p_a)/sqrt(length(p_a)),
                  sd(g_a)/sqrt(length(g_a)),
                  sd(p_g)/sqrt(length(p_g)))
  )
  rownames(out) <- c("P/A", "G/A", "P/G")

  out
}