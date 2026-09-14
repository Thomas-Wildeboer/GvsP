# Multivariate analysis of genetic and phenotypic variance-covariance matrices

Comparing simulated genetic and phenotypic variance-covariance matrices across a range of sample sizes, heritabilities, and pleiotropy levels.

## Workflow

1. Generate pedigree
2. Generate phenotypes along pedigree (`phensim.r`)
3. Estimate G matrix (WOMBAT)
4. Compare G, A, and P using utility scripts in `/util/`
5. Generate figures (`figures.Rmd`)

## Scripts

- **`main.r`** — Entry point. Sets simulation parameters (number of traits, sires, heritability, error structure) and calls `matrixsim()` for each parameter combination.
- **`matrixsim.r`** — Core simulation function. Generates a pedigree, simulates phenotypes via `phensim()`, estimates G via WOMBAT, then writes A (true breeding value covariance), G (WOMBAT-estimated), and P (phenotypic covariance) matrices to a CSV file. Results are saved as `*name*_matrices.csv` in the project root.
- **`phensim.r`** — Simulates phenotypes along a pedigree given A and E matrices. Based on the `phensim` function from the Pedantics package.
- **`figures.Rmd`** — Generates all ggplot figures from the matrices CSV files and analysis outputs.

## Data files

### Matrices CSVs (root directory)

- **`small_families_matrices.csv`** — Simulation results for small family structures: 20/60/200 families (10/30/100 sires × 2 dams/sire × 2 offspring/dam).
- **`large_families_matrices.csv`** — Simulation results for large family structures: 50/150/500 families (10/30/100 sires × 5 dams/sire × 3 offspring/dam).

Each row is one simulation replicate. Columns:
- `A` — True additive genetic covariance matrix (from phensim breeding values), stored as a space-separated vector of the lower triangle
- `G` — Estimated G matrix from WOMBAT, same format
- `P` — Phenotypic covariance matrix from phensim, same format
- `n` — Number of traits
- `h` — Heritability
- `multicoll` — Multicollinearity level of the true G matrix (`low`, `medium`, `high`)
- `f` — Number of families
- `e_type` — Environmental covariance structure
- `name` — Simulation run identifier

### Template files

- **`template_E.xlsx`**, **`template_n.xlsx`** — Excel templates used for parameter configurations.

## Directories

### `/util/`

Utility scripts that take the matrices CSV as input and compute comparison statistics, returning results as dataframes.

- `read_matrices.r` — Reads `*_matrices.csv` files into R
- `Angles.r` — Angles between matrix subspaces
- `AvgEvolvabilities.r` / `MaxEvolvabilities.r` — Average and maximum evolvability
- `Disparity.r` — Matrix disparity
- `KrzCorsFirst.r` / `KrzCorsLast.r` — Krzanowski subspace correlations (first and last axes)
- `MantelCors.r` / `MantelTests.r` — Mantel correlations and tests
- `RDifferences.r` — Response differences
- `RSkewers.r` — Random skewers analysis

### `/tables/`

CSV and XLSX outputs from the utility scripts

### `/figures/`

Intermediate PDF figures organised by family size:

- `small_families/` — Evolvabilities, Krzanowski correlations, Mantel correlations, random skewers, response differences, and regression plots for small family simulations
- `large_families/` — Same plots for large family simulations

### `/ms_figures/`

Final manuscript figures (PDFs)

### `/ms_tables/`

Final manuscript tables (XLSX)

### `/temp/`

Temporary files generated during WOMBAT runs.

## Dependencies

- R with packages: tidyverse, MCMCglmm, mixtools, monomvn, here, Matrix, evolvability, evolqg, reshape2
- [WOMBAT](http://didgeridoo.une.edu.au/km/wombat.php) available in the system path as `wombat`
- UNIX-like system (Linux, macOS, WSL on Windows)

## Acknowledgements

- Dr. J Sztepanacz (principal investigator)
- Alexander MacKenzie (code to interface with WOMBAT)
- Dr. M Videlier (helpful suggestions and discussions)
- Dr. M Morrissey (phensim function)
