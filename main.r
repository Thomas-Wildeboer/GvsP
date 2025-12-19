source("./matrixsim.r")

# set parameters
# 500 families: 100 sires, 5 dams per sire, 3 offspring per dam
# 150 families: 30 sires, 5 dams per sire, 3 offspring per dam
# 50 families: 10 sires, 5 dams per sire, 3 offspring per dam

# 200 families: 100 sires, 2 dams per sire, 2 offspring per dam
# 60 families: 30 sires, 2 dams per sire, 2 offspring per dam
# 20 families: 10 sires, 2 dams per sire, 2 offspring per dam

name <- "small_families_n8_f20_02"
iterations <- 25
n_traits <- c(8)
n_sires <- c(10)
dams_per_sire <- 2
offspring_per_dam <- 2
h <- c(0.1)
e_type <- c("ind")

# all combinations of parameters
params <- expand.grid(n_traits, n_sires, h, e_type)

# run simulations
for (i in 1:nrow(params)) {
    matrixsim(n_traits = params[i, 1], n_sires = params[i, 2], dams_per_sire = dams_per_sire, offspring_per_dam = offspring_per_dam, h = params[i, 3], e_type = params[i, 4], name = name, iterations = iterations)
}