# Generates pedigree, phenotypes, P matrix, then runs WOMBAT to estimate G matrix
# Finally, it writes the estimated G matrix from WOMBAT and the P matrix and true G matrix from R to csv files and deletes the WOMBAT output files
# Results are saved in the *name*_matrices.csv file
library(mixtools)
library(monomvn)
library(Matrix)
library(here)
source("./phensim.r")

# function that is required to fill the upper triangle of the matrix
upper.diag <- function(x) {
    m <- (-1+sqrt(1+8*length(x)))/2
    X <- lower.tri(matrix(NA, m, m), diag=TRUE)
    X[X==TRUE] <- x
    t(X)
}

# a set of different multicollinearity matrices are produced for each parameter combination
multicoll_levels <- c("low", "medium", "high")

# h -> heritability
# e_type -> environmental effects distribution reletaive to genetic effects (independent: "ind", similar/proportional: "sim", opposite: "opp")
# iterations -> number of iterations to run per parameter combination
# name -> name of the output file
`matrixsim` <-
function(n_traits=4, n_sires=5, dams_per_sire=3, offspring_per_dam=5, h=0.5, e_type="ind", iterations=1, name) {
    # set working directory
    dir.create(paste0(here(), "/temp/", name, "/"), showWarnings = FALSE)
    working_folder <- paste0(here(), "/temp/", name, "/")
    setwd(dir = working_folder)

    csv_filename <- paste0(working_folder, "../../", name , "_matrices.csv")

    # create csv file to store results if it doesn't exist
    if (!file.exists(csv_filename)) {
        df <- data.frame(matrix(ncol = 9, nrow = 0))
        x <- c("A", "G", "P", "n", "h", "multicoll", "f", "e_type", "name")
        colnames(df) <- x
        write.csv(df, csv_filename, row.names = FALSE)
    }
    
    # pedigree constants
    n_dams <- n_sires * dams_per_sire
    offspring_per_sire <- offspring_per_dam * dams_per_sire
    n_offspring <- n_dams * offspring_per_dam
    n_individuals <- n_sires + n_dams + n_offspring

    # create pedigree scaffold
    id_vector <- rep(NA, n_individuals) # IDs for all individuals
    no_parents <- rep(NA, n_sires + n_dams) # no parents for sires and dams in ped

    # set pedigree IDs
    for (i in 1:n_individuals) {
        id_vector[i] <- sprintf("%d", i)
    }

    # set offspring sire IDs
    offspring_sire_ids <- as.character(rep(1:n_sires, each=offspring_per_sire))

    # set offspring dam IDs
    offspring_dam_ids <- as.character(rep(n_sires+1:n_dams, each=offspring_per_dam))

    # sire and dam IDs for all individuals
    sire_vector <- c(no_parents, offspring_sire_ids)
    dam_vector <- c(no_parents, offspring_dam_ids)

    # ped for wombat (replace NA with 0)
    ped_wombat <- as.data.frame(cbind(id_vector, sire_vector, dam_vector))
    ped_wombat[is.na(ped_wombat)] <- 0

    # ped for phensim
    ped_phensim <- as.data.frame(cbind(id_vector, dam_vector, sire_vector))

    i <- 1

    while (i <= iterations) {
        print(paste0("Iteration ", i))

        # simulate additive effects and take covariance (this base matrix shared between pairs of high vs low multicollinearity)
        randomA <- cov(mvrnorm(n_individuals, mu=rep(0,n_traits), Sigma=diag(n_traits)))

        # create dataframes to store matrices
        low_df <- data.frame()
        medium_df <- data.frame()
        high_df <- data.frame()

        for (multicoll in multicoll_levels) {
            # modify eigenvalues of randomA to set multicollinearity
            # low multicollinearity: identity matrix
            # medium multicollinearity: geometrically decreasing coefficients with constant total sum
            # high: nearly all variance loaded on first principal component
            if (n_traits == 4) {
                if (multicoll == "high") {
                    eigenval_coeffs <- c(3.7, 0.1, 0.1, 0.1)
                } else if (multicoll == "medium") {
                    eigenval_coeffs <- c(2.2476931, 1.0432867, 0.4842508, 0.2247693)
                } else if (multicoll == "low") {
                    eigenval_coeffs <- c(1, 1, 1, 1)
                }
            } else if (n_traits == 8) {
                if (multicoll == "high") {
                    eigenval_coeffs <- c(7.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
                } else if (multicoll == "medium") {
                    eigenval_coeffs <- c(2.4164210, 1.7390636, 1.2515791, 0.9007436, 0.6482522, 0.4665378, 0.3357606, 0.2416421)
                } else if (multicoll == "low") {
                    eigenval_coeffs <- c(1, 1, 1, 1, 1, 1, 1, 1)
                }
            } else if (n_traits == 12) {
                if (multicoll == "high") {
                    eigenval_coeffs <- c(10.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
                } else if (multicoll == "medium") {
                    eigenval_coeffs <- c(2.4664950, 2.0006502, 1.6227890, 1.3162942, 1.0676868, 0.8660337, 0.7024666, 0.5697923, 0.4621761, 0.3748853, 0.3040810, 0.2466495)
                } else if (multicoll == "low") {
                    eigenval_coeffs <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
                }
            }

            # randomize order of eigenvalue coefficients
            # this prevents genetic variances from being larger than desired as R orders eigenvalues by size
            eigenval_coeffs <- sample(eigenval_coeffs)

            # modify A matrix to set multicollinearity - resize its eigenvalues with above coefficients
            modifiedA <- eigen(randomA)$vectors %*% (diag(eigenval_coeffs) * diag(eigen(randomA)$values)) %*% t(eigen(randomA)$vectors)

            # random environmental effects (E matrix, adjusted for heritability, independent, similar, or opposite to A matrix)
            if (e_type == "ind") {
                randomE <- diag(modifiedA) * ((1 - h) / h)
                randomE <- diag(randomE)
            } else if (e_type == "sim") {
                # adjust elements for heritability
                randomE <- modifiedA * ((1 - h) / h)
            } else if (e_type == "opp") {
                # adjust elements for heritability
                randomE <- modifiedA * ((1 - h) / h)

                # multiply lower diagonal elements by -1 to flip directions
                randomE[lower.tri(randomE, diag = FALSE)] <- randomE[lower.tri(randomE, diag = FALSE)] * -1

                # set upper diagonal elements to lower diagonal elements (to make symmetric)
                randomE[upper.tri(randomE, diag = FALSE)] <- randomE[lower.tri(randomE, diag = FALSE)]

                # ensure matrix is positive definite
                randomE <- nearPD(randomE, keepDiag = TRUE)$mat

                # convert to matrix type
                randomE <- as.matrix(randomE)
            }

            # simulate phenotypes -> dataframe
            phensim_out <- as.data.frame(phensim(ped_phensim, n_traits, modifiedA, randomE))

            # get phenotypes from phensim output
            phenotypes <- phensim_out[1:n_traits+1]

            # get breeding values from the rest of the phensim output
            breeding_values <- phensim_out[(n_traits+1):(n_traits*2)+1]

            # scale phenotypes and breeding values on same scale, so that V_P = 1 and V_A = h^2
            phenotypes <- scale(phenotypes, center = TRUE, scale = TRUE)
            breeding_values <- scale(breeding_values, center = TRUE, scale = TRUE) * sqrt(h)

            # P and true G (written here as A) matrices
            p_matrix <- cov(phenotypes)
            a_matrix <- cov(breeding_values)

            # print matrices for reference
            print("G matrix:")
            print(a_matrix)

            print("E matrix:")
            print(p_matrix-a_matrix)

            print("P matrix:")
            print(p_matrix)

            # create starting G for WOMBAT parameter file
            g_for_parameter <- cov(phenotypes)*h

            # making upper tri for starting G in parameter file
            g_for_parameter[upper.tri(g_for_parameter, diag = FALSE)] <- NA
            upper_tri_g <- as.data.frame(as.character(na.omit(as.vector(g_for_parameter))))

            # stacked data frame for WOMBAT
            stacked <- array(numeric(), c(nrow(ped_phensim)*n_traits, 3))

            # create stacked dataset for WOMBAT
            for (j in 1:nrow(ped_phensim)) {
                for (k in 1:n_traits) {
                    stacked[(j*n_traits)+k-(n_traits),] <- c(k, ped_phensim[j,1], phenotypes[j,k])
                }
            }

            # write stacked datafile and pedigree for WOMBAT
            write.table(stacked, paste0(working_folder, name, "Data.d"), sep = ' ', quote = FALSE, col.names = FALSE, row.names = FALSE)
            write.table(ped_wombat, paste0(working_folder, "ranPed.d"), sep = ' ', quote = FALSE, col.names = FALSE, row.names = FALSE)

            # shape parameter file
            parDat <- data.frame(nrow=1000)

            # first part of parameter file
            for (m in 1:n_traits){
                k <- (m-1)*3
                parDat[k+1,] <- paste0("tr", m, "  ", "traitno ", n_traits)
                parDat[k+2,] <- paste0("tr", m, "  ", "animal 0")
                parDat[k+3,] <- paste0("tr", m, "  ", "d", m, " 0")
            }

            # second part of parameter file
            tr <- data.frame(nrow=n_traits)

            # reduce rank model for 12 traits (ensures convergence)
            if (n_traits == 12) {
                reduction <- 1 # number of traits to reduce by

                # further reduce if h=0.1
                if (h == 0.1) {
                    reduction <- reduction + 1
                }

                # further reduce again if multicoll is high/medium
                if (multicoll == "high" || multicoll == "medium") {
                    reduction <- reduction + 1
                }

                unchanging_variables_par <-
                    c(paste0("ANAL MUV PC ", n_traits), # "PC" after "ANAL MUV" specifies reduced rank model
                    c("PEDS ranPed.d"),
                    c("END DATA"),
                    c("MODEL"),
                    c("RAN animal NRM"),
                    c("END MODEL"))
                
                # creating dynamic variables for parameter file
                for (t in 1:n_traits){
                    tr[t,] <- paste0("tr d", t, " ", t)
                    var_animal <- paste0("VAR animal ", t, " ", t-reduction) # reduce second t for reduced rank model
                    var_resid <- paste0("VAR resid ", t, " ", t)
                }
            } else {
                unchanging_variables_par <-
                    c(paste0("ANAL MUV ", n_traits), # full rank model
                    c("PEDS ranPed.d"),
                    c("END DATA"),
                    c("MODEL"),
                    c("RAN animal NRM"),
                    c("END MODEL"))

                # creating dynamic variables for parameter file
                for (t in 1:n_traits){
                    tr[t,] <- paste0("tr d", t, " ", t)
                    var_animal <- paste0("VAR animal ", t, " ", t)
                    var_resid <- paste0("VAR resid ", t, " ", t)
                }
            }

            # create stacked parameter file and output
            out <- unname(rbind(unchanging_variables_par[1], unchanging_variables_par[2], paste0("DATA ", name, "Data.d"), parDat, unchanging_variables_par[3], unchanging_variables_par[4], unchanging_variables_par[5], tr, unchanging_variables_par[6], var_animal))
            out2 <- rbind(upper_tri_g, var_resid, upper_tri_g)

            colnames(out) <- colnames(out2)

            out3 <- rbind(out, out2)

            write.table(out3, paste0(name, ".par"), row.names = FALSE, quote = FALSE, col.names = paste0("###G", i, " ", "dim=", n_traits)) 

            # run wombat assuming it is installed with "wombat" command in path
            system(paste0("wombat --cycle --pxai ", name, ".par"))

            # read iterates file
            iterates <- paste(readLines(paste0(working_folder, "Iterates")), collapse="\n")

            # check if "converged" is in iterates -> if not, need to try again
            if (!grepl("converged", iterates)) {
                print("WOMBAT did not converge")

                # delete WOMBAT files and exit loop
                system(paste0("rm -r ", working_folder, "*"))

                break
            }

            # grab G estimate from WOMBAT output
            bestpoint <- scan(paste0(working_folder, "BestPoint"), skip=1)

            matrix_size <- ((n_traits*(n_traits+1))/2)+1 # calculate the size of the matrix+1 (+1 needed to start at correct position)

            g_matrix <- upper.diag(bestpoint[matrix_size:length(bestpoint)]) # additive cov matrix is after the residual cov matrix
            g_matrix <- Matrix::forceSymmetric(g_matrix) # make into symmetric matrix
            g_matrix <- as.matrix(g_matrix) # convert to matrix

            # delete WOMBAT files
            system(paste0("rm -r ", working_folder, "*"))

            # add matrix to df object
            if (multicoll == "low") {
                low_df <- data.frame(paste(unlist(a_matrix), collapse=" "), paste(unlist(g_matrix), collapse=" "), paste(unlist(p_matrix), collapse=" "), n_traits, h, multicoll, n_dams, e_type, name)
            } else if (multicoll == "medium") {
                medium_df <- data.frame(paste(unlist(a_matrix), collapse=" "), paste(unlist(g_matrix), collapse=" "), paste(unlist(p_matrix), collapse=" "), n_traits, h, multicoll, n_dams, e_type, name)            
            } else if (multicoll == "high") {
                high_df <- data.frame(paste(unlist(a_matrix), collapse=" "), paste(unlist(g_matrix), collapse=" "), paste(unlist(p_matrix), collapse=" "), n_traits, h, multicoll, n_dams, e_type, name)
            }
        }

        # write matrices to csv if matrices from all levels of multicollinearity were successfully estimated
        if (nrow(low_df) > 0 && nrow(medium_df) > 0 && nrow(high_df) > 0) {
            write.table(low_df, file = paste0(here(), "/", name, "_matrices.csv"), sep = ",", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)
            write.table(medium_df, file = paste0(here(), "/", name, "_matrices.csv"), sep = ",", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)
            write.table(high_df, file = paste0(here(), "/", name, "_matrices.csv"), sep = ",", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)
        } else {
           # clear df objects
            rm(low_df)
            rm(medium_df)
            rm(high_df)

            # retry iteration if one of the matrix estimations failed to converge
            next
        }

        # clear df objects
        rm(low_df)
        rm(medium_df)
        rm(high_df)

        i <- i+1 # increment iteration
    }
}
