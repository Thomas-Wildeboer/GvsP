MantelTests <- function(a_list, g_list, p_list, alpha = 0.05) {
    # helper to extract p-values from MantelCor output
    get_p <- function(res) res[2]

    # A vs P
    mantel_p <- Map(function(p, a) MantelCor(cov2cor(p), cov2cor(a)),
                    p_list, a_list)
    pvals_p <- unlist(lapply(mantel_p, get_p))

    # A vs G
    mantel_g <- Map(function(g, a) MantelCor(cov2cor(g), cov2cor(a)),
                    g_list, a_list)
    pvals_g <- unlist(lapply(mantel_g, get_p))

    # P vs G
    mantel_pg <- Map(function(g, p) MantelCor(cov2cor(g), cov2cor(p)),
                    g_list, p_list)
    pvals_pg <- unlist(lapply(mantel_pg, get_p))

    # compute proportions and SEs
    N <- length(pvals_p)
    
    prop_p  <- mean(pvals_p  < alpha)
    prop_g  <- mean(pvals_g  < alpha)
    prop_pg <- mean(pvals_pg < alpha)

    se_p  <- sqrt(prop_p  * (1 - prop_p)  / N)
    se_g  <- sqrt(prop_g  * (1 - prop_g)  / N)
    se_pg <- sqrt(prop_pg * (1 - prop_pg) / N)

    # return base R data.frame
    df <- data.frame(
        Proportion = c(prop_p, prop_g, prop_pg),
        SE         = c(se_p,  se_g,  se_pg),
        row.names  = c("P,A", "G,A", "P,G"),
        stringsAsFactors = FALSE
    )

    return(df)
}