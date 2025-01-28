source("./matrixsim.r")

# set parameters
# 500 families: 100 sires, 5 dams per sire, 3 offspring per dam
# 150 families: 30 sires, 5 dams per sire, 3 offspring per dam
# 50 families: 10 sires, 5 dams per sire, 3 offspring per dam

name <- "example"
iterations <- 100
n_traits <- c(4, 8, 12)
n_sires <- c(10, 30, 100)
dams_per_sire <- 5
offspring_per_dam <- 3
h <- c(0.1, 0.5)
e_type <- c("ind", "sim", "opp")

# all combinations of parameters
params <- expand.grid(n_traits, n_sires, h, e_type)

# run simulations
for (i in 1:nrow(params)) {
    matrixsim(n_traits = params[i, 1], n_sires = params[i, 2], dams_per_sire = dams_per_sire, offspring_per_dam = offspring_per_dam, h = params[i, 3], e_type = params[i, 4], name = name, iterations = iterations)
}