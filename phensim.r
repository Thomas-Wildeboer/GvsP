# based on phensim function from Pedantics R/CRAN package by Michael Morrissey <michael.morrissey@st-andrews.ac.uk>
# generates phenotypes with given pedigree and provided random A/E matrices

library(MCMCglmm)

`phensim` <-
function(pedigree, traits = 1, randomA = NULL, randomE = NULL) {
    # Simulate breeding values and environmental effects
    n <- length(pedigree[,1])

    breedingValues <- rbv(pedigree, randomA) # rbv gives breeding values
    
    environEffects <- rmvnorm(n, rep(0, dim(randomE)[2]), randomE) # random environmental effects from multivariate normal distribution

    effectsTable <- as.data.frame(cbind(
        pedigree, breedingValues, environEffects
    ))

    for (x in 1:traits) {
        names(effectsTable)[3+x] <- paste("a_tr", x, sep="")
    }

    for(x in 1:traits) {
        names(effectsTable)[3+dim(randomA)[2]+x] <- paste("e_tr", x, sep="")
    }

    calcPhen <- effectsTable[, 4:length(effectsTable[1,])]
    retainIndex<-c(rep(1, traits), rep(1, traits))
    retainIndex <- matrix(retainIndex, length(calcPhen[,1]), length(calcPhen[1,]), byrow=TRUE)
    calcPhen <- calcPhen * retainIndex

    for (x in 1:traits) {
        ind <- array(0, dim=traits)
        ind[x] <- 1
        ind <- matrix(ind, length(calcPhen[,1]), length(calcPhen[1,]), byrow=TRUE)
        calcPhenTrait <- calcPhen * ind
        effectsTable[, (dim(effectsTable)[2]+1)] <- rowSums(calcPhenTrait)
        names(effectsTable)[dim(effectsTable)[2]] <- paste("Phen_tr", x, sep="")
    }

    phenotypes <- as.data.frame(effectsTable[,(dim(effectsTable)[2]-traits+1):dim(effectsTable)[2]])
    phenotypes <- cbind(pedigree[,1], phenotypes)

    names(phenotypes)[1] <- "id"

    for (x in 1:traits) {
        names(phenotypes)[x+1] <- paste("trait_", x, sep="")
    }

    # Output results
    output <- list(phenotypes=phenotypes, breedingValues=breedingValues)
    output
}