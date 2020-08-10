# Use the pedigreemm library to run a pedigree BLUP
library(AlphaSimR)
library(pedigreemm)

#Simulation
founderPop = runMacs(100,8,100)

SP = SimParam$new(founderPop)
SP$addTraitAD(100, meanDD=0.2, varDD=0.1)
SP$setVarE(h2=0.2)
SP$setTrackPed(TRUE)
SP$setGender("yes_sys")

trainPop = pop = newPop(founderPop)
for(i in 1:10){
  pop = selectCross(pop,nMale=5,nFemale=40,nCrosses=100)
  trainPop = c(trainPop, pop)
}

#Construct pedigree
id = as.factor(1:SP$lastId)
dam = SP$pedigree[,1]
dam[dam==0L] = NA
sire = SP$pedigree[,2]
sire[sire==0L] = NA
ped = pedigree(dam,sire,id)

#Run model
df = data.frame(y=pheno(pop),ID=id[as.numeric(pop@id)])
ans = pedigreemm(y~(1|ID),data=df,
                 pedigree=list(ID=ped),
                 control=lmerControl(
                   check.nobs.vs.nlev="ignore",
                   check.nobs.vs.nRE="ignore"))

summary(ans)

cor(bv(pop), bv(trainPop)[1001:1100])
cor(ranef(ans)$ID, bv(trainPop)[1001:1100])
cor(ranef(ans)$ID, bv(pop))


