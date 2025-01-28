# Mantel tests of correlations between matrices

MantelTests <- function(a_list, g_list, p_list) {
    # A vs P ------------------------------------------------------------------
    mantel_p <- map2(p_list, a_list, ~MantelCor(cov2cor(.x), cov2cor(.y))) # mantel cor requires correlation matrices

    pvals_p <- lapply(mantel_p, function(x) x[2])
    pvals_p <- as.numeric(pvals_p)

    # A vs G ------------------------------------------------------------------
    mantel_g <- map2(g_list, a_list, ~MantelCor(cov2cor(.x), cov2cor(.y)))

    pvals_g <- lapply(mantel_g, function(x) x[2])
    pvals_g <- as.numeric(pvals_g)
    
    # P vs G ------------------------------------------------------------------
    mantel_pg <- map2(g_list, p_list, ~MantelCor(cov2cor(.x), cov2cor(.y)))
    
    pvals_pg <- lapply(mantel_pg, function(x) x[2])
    pvals_pg <- as.numeric(pvals_pg)

    # create data frame of results
    pvalues <- data.frame(
        c(mean(pvals_p), mean(pvals_g), mean(pvals_pg))
    )

    # set column names
    colnames(pvalues) <- c(
        "Mean Mantel P-Value"
    )

    # set row names
    rownames(pvalues) <- c(
        "P,A",
        "G,A",
        "P,G"
    )

    return(pvalues)
}
