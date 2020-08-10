# This script illustrates heterosis and inbreeding depression
library(AlphaSimR)

# Create founder haplotypes for two breeds
# The split parameter models the number generations ago that the breeds split.
# This will take some time.
founderPop = runMacs(nInd = 1000,
                     nChr = 10,
                     segSites = 1000,
                     split = 100) 

# Set simulation parameters
SP = SimParam$new(founderPop)

# Add a trait with dominance
SP$addTraitAD(nQtlPerChr = 1000,
              meanDD = 0.2, # Mean value of dominance degrees 
              varDD = 0.1) # Variance of dominance degrees

# Create a populations for each breed
# Note that AlphaSimR models the breeds as equal size
# This means individuals 1 to 500 belong to one breed and
# individuals 501 to 1000 belong to another breed
breedA = newPop(founderPop[1:500])
breedB = newPop(founderPop[501:1000])

# Create fully inbred populations from each breed using makeDH
inbredA = makeDH(breedA)
inbredB = makeDH(breedB)

# Estimate inbreeding depression for each breed using the inbreds
# Note that this value gets larger with more QTL and/or higher meanDD
meanG(breedA) - meanG(inbredA)
meanG(breedB) - meanG(inbredB)

# Randomly create 1000 crossbred animals
crossbred = randCross2(breedA, breedB, nCrosses = 1000)

# Estimate heterosis using the crossbred animals
# Note that this value gets larger with more QTL, higher meanDD, 
# and/or higher value for split in the runMacs command
meanG(crossbred) - (meanG(breedA) + meanG(breedB))/2


