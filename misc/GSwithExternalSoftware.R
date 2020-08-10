library(AlphaSimR)
library(rrBLUP)

founderPop = runMacs(100,10,300) #100 animals with 10 chromosomes and 300 sites per chromosome

SP = SimParam$new(founderPop)
SP$addTraitA(100) #Add a trait with 100 QTL per chromosome
SP$setVarE(h2=0.5) #Set trait heritability to 0.5
SP$addSnpChip(200) #Add a SNP chip with 200 markers per chromosome

pop = newPop(founderPop) #Create initial population

y = pheno(pop) #Access phenotype data
M = pullSnpGeno(pop) #Access SNP genotype data
M = M-1 #Convert markers from 0,1,2 coding to -1,0,1 coding
G = A.mat(M) #Create a G using an rrBLUP function
ans = mixed.solve(y=y, K=G) #Run a GBLUP using the rrBLUP solver
EBV = ans$u #Extract EBVs from rrBLUP output

# Examine prediction accuracy
cor(EBV, gv(pop))

# Insert externally created EBV into the population
pop@ebv = as.matrix(EBV) #Assign EBVs to population 
