# Mantel correlations between matrices

MantelCors <- function(a_list, g_list, p_list) {
    # A vs P ------------------------------------------------------------------
    mantel_p <- map2(p_list, a_list, ~MantelCor(cov2cor(.x), cov2cor(.y))) # mantel cor requires correlation matrices

    cors_p <- lapply(mantel_p, function(x) x[1])
    cors_p <- as.numeric(cors_p)

    cors_p_se <- sd(cors_p)/sqrt(length(cors_p))

    # A vs G ------------------------------------------------------------------
    mantel_g <- map2(g_list, a_list, ~MantelCor(cov2cor(.x), cov2cor(.y)))

    cors_g <- lapply(mantel_g, function(x) x[1])
    cors_g <- as.numeric(cors_g)

    cors_g_se <- sd(cors_g)/sqrt(length(cors_g))
    
    # P vs G ------------------------------------------------------------------
    mantel_pg <- map2(g_list, p_list, ~MantelCor(cov2cor(.x), cov2cor(.y)))
    
    cors_pg <- lapply(mantel_pg, function(x) x[1])
    cors_pg <- as.numeric(cors_pg)
    
    cors_pg_se <- sd(cors_pg)/sqrt(length(cors_pg))

    # create data frame of results
    mantelcors <- data.frame(
        c(mean(cors_p), mean(cors_g), mean(cors_pg)),
        c(cors_p_se, cors_g_se, cors_pg_se)
    )

    # set column names
    colnames(mantelcors) <- c(
        "Mean Mantel Correlation",
        "SE"
    )

    # set row names
    rownames(mantelcors) <- c(
        "P,A",
        "G,A",
        "P,G"
    )

    return(mantelcors)
}
