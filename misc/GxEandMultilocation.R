# This script generates multi-environment data with GxE

library(AlphaSimR)

founderPop = runMacs(100, 10, 1000)

SP = SimParam$
  new(founderPop)$
  addTraitAG(1000,varGxE=2)

pop = newPop(founderPop)

# Measure location 1 phenotypes
p = runif(1) # A p-value for the environmnetal covariate
loc1 = setPheno(pop, varE=4, p=p, onlyPheno=TRUE)

# Measure location 2 phenotypes
p = runif(1) # A p-value for the environmnetal covariate
loc2 = setPheno(pop, varE=4, p=p, onlyPheno=TRUE)

# Measure location 3 phenotypes
p = runif(1) # A p-value for the environmnetal covariate
loc3 = setPheno(pop, varE=4, p=p, onlyPheno=TRUE)

# Calculate across location average
allLoc = cbind(loc1,loc2,loc3)
allLocMean = rowMeans(allLoc)

# Add average phenotype to population 
# Note that phenotypes must be matrices
# The rows of the matrix equal the individuals
# and the columns are the traits
pop@pheno = as.matrix(allLocMean)

# Select the best 10 individuals
selPop = selectInd(pop, 10)
meanG(selPop)




