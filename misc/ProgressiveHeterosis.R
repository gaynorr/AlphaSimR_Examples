# This script illustrates progressive heterosis in polyploids

library(AlphaSimR)

het2sc = het4sc = het2dc = het4dc = numeric(100)
for(REP in 1:100){
  founderPop = quickHaplo(nInd=4000, nChr=4, segSites=100, ploidy=4)
  
  SP = SimParam$
    new(founderPop)$
    addTraitAD(100, meanDD=0.6, varDD=0.2)$
    setVarE(h2=0.3)
  
  # Create 4 tetraploid populations
  tetraploid = vector("list", 4)
  tetraploid[[1]] = newPop(founderPop[1:100])
  tetraploid[[2]] = newPop(founderPop[101:200])
  tetraploid[[3]] = newPop(founderPop[201:300])
  tetraploid[[4]] = newPop(founderPop[301:400])
  
  # Create diploid populations from tetraploid populations
  diploid = vector("list", 4)
  for(i in 1:4){
    diploid[[i]] = reduceGenome(tetraploid[[i]])
  }
  
  # Perform 5 generations of recurrent selection
  for(gen in 1:5){ 
    for(i in 1:4){
      tetraploid[[i]] = selectCross(tetraploid[[i]], nInd=10, nCrosses=100)
      diploid[[i]] = selectCross(diploid[[i]], nInd=10, nCrosses=100)
    }
  }
  
  # Create single cross populations (1/2 and 3/4)
  tetraSC = list(randCross2(tetraploid[[1]], tetraploid[[2]], nCrosses=100),
                 randCross2(tetraploid[[3]], tetraploid[[4]], nCrosses=100))
  diploidSC = list(randCross2(diploid[[1]], diploid[[2]], nCrosses=100),
                   randCross2(diploid[[3]], diploid[[4]], nCrosses=100))
  
  # Create double cross populations
  tetraDouble = randCross2(tetraSC[[1]], tetraSC[[2]], nCrosses=1000)
  diploidDouble = randCross2(diploidSC[[1]], diploidSC[[2]], nCrosses=1000)
  
  
  # Calculate mean of all base populations
  meanTetra = mean(unlist(lapply(tetraploid,meanG)))
  meanDiploid = mean(unlist(lapply(diploid,meanG)))
  
  # Calculate mean of single cross populations
  meanTetraSC = mean(unlist(lapply(tetraSC,meanG)))
  meanDiploidSC = mean(unlist(lapply(diploidSC,meanG)))
  
  # Measure single cross heterosis
  het4sc[REP] = meanTetraSC - meanTetra
  het2sc[REP] = meanDiploidSC - meanDiploid
  
  # Measure double cross heterois
  het4dc[REP] = meanG(tetraDouble) - meanTetra 
  het2dc[REP] = meanG(diploidDouble) - meanDiploid 
}
boxplot(het2dc,het2sc,het4dc,het4sc)
boxplot(het2dc-het2sc,het4dc-het4sc)

