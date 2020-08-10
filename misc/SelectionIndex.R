# A demonstration of using a selection index

## Initial parameters
# Genetic correlation matrix
genCor = matrix(c( 1,   -0.6,  0.1,
                  -0.6,  1,    0.3,
                   0.1,  0.3,  1), ncol=3)

# Trait heritabilities
h2 = c(0.3, 0.5, 0.3)

# Trait economic weights
w = c(3, -1, 1)

## Simulation using above parameters
library(AlphaSimR)

founderPop = quickHaplo(nInd=1000, nChr=10, segSites=100)

SP = SimParam$
  new(founderPop)$
  addTraitA(100, mean=c(0,0,0), var=c(1,1,1),
            corA=genCor)$
  setVarE(h2=h2)

pop = newPop(founderPop)

# Calculation average true economic value in the original population
mean(gv(pop)%*%w)

## Naive selection index, apply economic weights to phenotypes
naive = selectInd(pop, nInd=100, 
                  trait=selIndex, # Index function
                  b=w) # Define weights as economic weights

# Calculation average true economic value in the naive population
mean(gv(naive)%*%w)


## Smith-Hazel selection index
# Calculate correct weights assuming variances are known
b = smithHazel(econWt=w, varG=varG(pop), varP=varP(pop))
best = selectInd(pop, nInd=100, 
                  trait=selIndex, # Index function
                  b=b) # Define weights as Smith-Hazel weights

# Calculation average true economic value in the best population
mean(gv(best)%*%w) 


