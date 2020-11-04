# Performance test script
# Measures timing and memory use (Windows only)
# for a simulation modeling 100 generations of 
# recurrent selection. Each generation consists 
# of 1000 parents and 10,000 offspring. Each 
# individual has 10 chromosomes with 10,000
# loci per chromosome.
# NOTE: memory.size only works on Windows

library(AlphaSimR)

system.time({
  founderPop = quickHaplo(nInd = 1000,
                          nChr = 10,
                          segSites = 1000)
  
  SP = SimParam$new(founderPop)$
    addTraitA(nQtlPerChr = 1000)$
    setVarE(h2 = 0.3)
  
  pop = newPop(founderPop)
  
  for(generation in 1:100){
    pop = selectCross(pop = pop, 
                      nInd = 1000,
                      nCrosses = 10000)
  }
})

memory.size(TRUE)

