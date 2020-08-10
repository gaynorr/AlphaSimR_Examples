library(AlphaSimR)

founderPop = quickHaplo(nInd=1000, nChr=10, segSites=600)

SP = SimParam$
  new(founderPop)$
  addTraitA(100)$
  addSnpChip(500)$
  setVarE(h2=0.4)
  
# Split the founders into two population
pop1 = newPop(founderPop[1:500])
pop2 = newPop(founderPop[501:1000])

# Create a third population derived from the first
pop3 = randCross(pop1, nCrosses=1000)

# Train a GS model using the first population
gsModel = RRBLUP(pop1)

# Set EBVs for all populations using the GS model
pop1 = setEBV(pop1, gsModel)
pop2 = setEBV(pop2, gsModel)
pop3 = setEBV(pop3, gsModel)

# Measure prediction accuracy
cor(gv(pop1), ebv(pop1)) 
cor(gv(pop2), ebv(pop2)) 
cor(gv(pop3), ebv(pop3))

# Explainations for observed accuraccies
# 
# pop1: Accuracy high due to being the training population.
# pop2: Accuracy near zero due to no LD between SNP and QTL.
# There would be some accuracy if runMacs is used instead, 
# because MaCS will generate some LD between SNP and QTL.
# pop3: Reasonable accuracy due to LD between SNP and QTL. 
# The LD is due to the individuals being progeny of the 
# training population.
  


