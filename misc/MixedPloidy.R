# Simulates 1 breeding cycle for improving triploids 
# The triploids are created by crossing diploids to tetraploids
# Selection is applied to diploids and tetraploids using GCA
# GCA is calculated by evaluating the triploid progeny

library(AlphaSimR)

tetraFounder = runMacs(nInd=200, nChr=10, 
                       segSites=1000, ploidy=4)

SP = SimParam$
  new(tetraFounder)$
  addTraitAD(1000, meanDD=0.2, varDD=0.1)$
  setVarE(h2=0.3)

# Create tetraploid parents from first 100 founders
tetraploids = newPop(tetraFounder[1:100])

# Create diploid parents from the other founders
diploids = reduceGenome(newPop(tetraFounder[101:200]))

# Create triploid plants for evaluation
# You will get 5 crosses per parents, because balance=TRUE
triploids = randCross2(tetraploids, diploids, nCrosses=500)

# Analyze triploid data to calculate GCA for parents
gca = calcGCA(triploids)

# Match tetraploid ids to their GCA
tetraMatch = match(tetraploids@id, 
                   gca$GCAf$id) # GCA for females (tetraploids)
tetraGCA = gca$GCAf$Trait1[tetraMatch]

# Select the best 20 tetraploids based on their GCA
tetraploids2 = tetraploids[order(tetraGCA, decreasing=TRUE)[1:20]]

# Create 100 new tetraploids by crossing selected tetraploids
tetraploids2 = randCross(tetraploids2, 100)

# Match diploid ids to their GCA
diploidMatch = match(diploids@id, 
                   gca$GCAm$id) # GCA for males (diploids)
diploidGCA = gca$GCAm$Trait1[diploidMatch]

# Select the best 10 diploids based on their GCA
diploids2 = diploids[order(diploidGCA, decreasing=TRUE)[1:20]]

# Create 100 new diploids by crossing selected diploids
diploids2 = randCross(diploids2, 100)

# Create new triploids from the new diploids and tetraploids
triploids2 = randCross2(tetraploids2, diploids2, nCrosses=500)

# Measure the genetic gain for all ploidy levels
meanG(tetraploids2) - meanG(tetraploids)
meanG(diploids2) - meanG(diploids)
meanG(triploids2) - meanG(triploids)
