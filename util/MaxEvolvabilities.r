# Maximum eigenvalue of a covariance matrix
# Max evolvability over random selection gradients = max eigenvalue of G (Hansen and Houle 2008)
# Assuming selection gradients are uniformly distributed in n-dimensional space

library(purrr)

MaxEvolvabilities <- function(a_list, g_list, p_list) {
    # get average eigenvalues of G matrices
    g_eigenvalues <- map(g_list, ~eigen(.x)$values)
    g_eigenvalues <- map_dbl(g_eigenvalues, max)

    # get average eigenvalues of P matrices
    p_eigenvalues <- map(p_list, ~eigen(.x)$values)
    p_eigenvalues <- map_dbl(p_eigenvalues, max)

    # get average eigenvalues of A matrices
    a_eigenvalues <- map(a_list, ~eigen(.x)$values)
    a_eigenvalues <- map_dbl(a_eigenvalues, max)

    # relative max evolvabilities
    evolvs_p <- p_eigenvalues/a_eigenvalues
    evolvs_g <- g_eigenvalues/a_eigenvalues
    evolvs_pg <- p_eigenvalues/g_eigenvalues

    # get standard errors
    evolvs_p_se <- sd(evolvs_p)/sqrt(length(evolvs_p))
    evolvs_g_se <- sd(evolvs_g)/sqrt(length(evolvs_g))
    evolvs_pg_se <- sd(evolvs_pg)/sqrt(length(evolvs_pg))
    
    # create data frame of results
    evolvs <- data.frame(
        c(mean(evolvs_p), mean(evolvs_g), mean(evolvs_pg)),
        c(evolvs_p_se, evolvs_g_se, evolvs_pg_se)
    )

    # set column names
    colnames(evolvs) <- c(
        "Relative Evolvability",
        "SE"
    )

    # set row names
    rownames(evolvs) <- c(
        "P/A",
        "G/A",
        "P/G"
    )

    return(evolvs)
}
