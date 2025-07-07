# Response Distance from Hansen and Houle (2008)
library(evolvability)
library(purrr)

ResponseDifference <- function(A, X, i=1000) {
    # get random selection gradients
    betas <- randomBeta(n=i, k=dim(X)[1])

    # get differences between response vectors from A and X
    differences <- map_dbl(1:i, ~((norm((X-A) %*% betas[, .x]))/
                            norm(A %*% betas[, .x])))

    return(mean(differences))
}

RDifferences <- function(a_list, g_list, p_list) {
    # A vs P ------------------------------------------------------------------
    rdiff_p <- map2_dbl(a_list, p_list, ~ResponseDifference(.x, .y))
    rdiff_p_se <- sd(rdiff_p)/sqrt(length(rdiff_p))
    
    # A vs G ------------------------------------------------------------------
    rdiff_g <- map2_dbl(a_list, g_list, ~ResponseDifference(.x, .y))
    rdiff_g_se <- sd(rdiff_g)/sqrt(length(rdiff_g))

    # P vs G ------------------------------------------------------------------
    rdiff_gp <- map2_dbl(p_list, g_list, ~ResponseDifference(.x, .y))
    rdiff_gp_se <- sd(rdiff_gp)/sqrt(length(rdiff_gp))
    
    # create data frame of results
    rdiffs <- data.frame(
        c(mean(rdiff_p), mean(rdiff_g), mean(rdiff_gp)),
        c(rdiff_p_se, rdiff_g_se, rdiff_gp_se)
    )

    # set column names
    colnames(rdiffs) <- c(
        "Mean Response Difference",
        "SE"
    )

    # set row names
    rownames(rdiffs) <- c(
        "P,A",
        "G,A",
        "G,P"
    )

    return(rdiffs)
}