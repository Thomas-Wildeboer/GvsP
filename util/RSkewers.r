# Random skewer correlations between matrices

RSkewers <- function(a_list, g_list, p_list) {
    # A vs P ------------------------------------------------------------------
    rskewers_p <- map2(p_list, a_list, ~RandomSkewers(.x, .y, nperm = 1000))
    rskewers_p <- unlist(rskewers_p)
    rskewers_p <- as.data.frame(rskewers_p)
    
    # get correlations: only keep first of every 3 rows of rskewers
    rskewers_p_cor <- rskewers_p[seq(1, nrow(rskewers_p), 3), ]

    rskewers_p_cor_se <- sd(rskewers_p_cor)/sqrt(length(rskewers_p_cor))
    
    # A vs G ------------------------------------------------------------------
    rskewers_g <- map2(g_list, a_list, ~RandomSkewers(.x, .y, nperm = 1000))
    rskewers_g <- unlist(rskewers_g)
    rskewers_g <- as.data.frame(rskewers_g)
    
    # get correlations: only keep first of every 3 rows of rskewers
    rskewers_g_cor <- rskewers_g[seq(1, nrow(rskewers_g), 3), ]
    
    rskewers_g_cor_se <- sd(rskewers_g_cor)/sqrt(length(rskewers_g_cor))

    # P vs G ------------------------------------------------------------------
    rskewers_pg <- map2(g_list, p_list, ~RandomSkewers(.x, .y, nperm = 1000))
    rskewers_pg <- unlist(rskewers_pg)
    rskewers_pg <- as.data.frame(rskewers_pg)
    
    # get correlations: only keep first of every 3 rows of rskewers
    rskewers_gp_cor <- rskewers_pg[seq(1, nrow(rskewers_pg), 3), ]
    
    rskewers_gp_cor_se <- sd(rskewers_gp_cor)/sqrt(length(rskewers_gp_cor))
    
    # create data frame of results
    rskewers <- data.frame(
        c(mean(rskewers_p_cor), mean(rskewers_g_cor), mean(rskewers_gp_cor)),
        c(rskewers_p_cor_se, rskewers_g_cor_se, rskewers_gp_cor_se)
    )

    # set column names
    colnames(rskewers) <- c(
        "Mean RS Correlation",
        "SE"
    )

    # set row names
    rownames(rskewers) <- c(
        "P,A",
        "G,A",
        "G,P"
    )

    return(rskewers)
}
