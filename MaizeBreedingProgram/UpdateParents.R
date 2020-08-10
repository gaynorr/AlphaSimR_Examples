#Sets new parents for inbreds
#Uses YT4 and YT3 inbreds with best YT2 inbreds
MaleParents = c(MaleInbredYT4,MaleInbredYT3)
nSelYT2 = nParents - MaleParents@nInd
if(nSelYT2>0){
  MaleParents = c(MaleParents,selectInd(MaleYT2,nSelYT2))
}else if(nSelYT2<0){
  MaleParents = MaleParents[1:nParents]
}

FemaleParents = c(FemaleInbredYT4,FemaleInbredYT3)
nSelYT2 = nParents - FemaleParents@nInd
if(nSelYT2>0){
  FemaleParents = c(FemaleParents,selectInd(FemaleYT2,nSelYT2))
}else if(nSelYT2<0){
  FemaleParents = FemaleParents[1:nParents]
}
