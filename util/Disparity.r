# Disparity between correlation matrices (Willis 1991)
# Average of absolute differences between pairwise correlations as a percentage of the mean absolute genetic correlation

library(purrr)

# Takes correlation matrices
Disparity <- function(a_list, g_list, p_list) {
    p_a_disparity <- map2_dbl(p_list, a_list, ~{
        mean_gen_corr <- mean(abs(.y[lower.tri(.y)]))
        mean(abs(.x[lower.tri(.x)] - .y[lower.tri(.y)])) / mean_gen_corr * 100
    })
    
    g_a_disparity <- map2_dbl(g_list, a_list, ~{
        mean_gen_corr <- mean(abs(.y[lower.tri(.y)]))
        mean(abs(.x[lower.tri(.x)] - .y[lower.tri(.y)])) / mean_gen_corr * 100
    })
    
    p_g_disparity <- map2_dbl(p_list, g_list, ~{
        mean_gen_corr <- mean(abs(.y[lower.tri(.y)]))
        mean(abs(.x[lower.tri(.x)] - .y[lower.tri(.y)])) / mean_gen_corr * 100
    })

    # Standard errors
    g_a_disparity_se <- sd(g_a_disparity) / sqrt(length(g_a_disparity))
    p_a_disparity_se <- sd(p_a_disparity) / sqrt(length(p_a_disparity))
    p_g_disparity_se <- sd(p_g_disparity) / sqrt(length(p_g_disparity))

    # Create data frame of results
    disparities <- data.frame(
        Disparity = c(mean(p_a_disparity), mean(g_a_disparity), mean(p_g_disparity)),
        SE = c(p_a_disparity_se, g_a_disparity_se, p_g_disparity_se)
    )

    # Set row names
    rownames(disparities) <- c("P/A", "G/A", "P/G")

    return(disparities)
}
