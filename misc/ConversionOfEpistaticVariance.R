# An example showing epistatic variance being "converted" to 
# additive genetic variance

library(AlphaSimR)
library(ggplot2)

## User set variances
# VarA = 1, always
VarD = 1 
VarAA = 1
VarE = 1
n = 1000 # Population size, must be divisable by 10
# Note that conversion is due to drift, so it occurs
# sooner in smaller populations

## Simultation
founderPop = quickHaplo(nInd=n,nChr=10,segSites=100)

SP = SimParam$
  new(founderPop)$
  addTraitADE(100, varDD=VarD*2, relAA=VarAA)$
  setVarE(varE=VarE)

pop = newPop(founderPop)
MEAN = VARG = VARA = VARD = VARAA = numeric(101)
tmp = genParam(pop)
MEAN[1] = tmp$mu
VARG[1] = tmp$varG
VARA[1] = tmp$varA
VARD[1] = tmp$varD
VARAA[1] = tmp$varAA
for(i in 2:101){
  pop = selectCross(pop,n*0.1,nCrosses=n)
  tmp = genParam(pop)
  MEAN[i] = tmp$mu
  VARG[i] = tmp$varG
  VARA[i] = tmp$varA
  VARD[i] = tmp$varD
  VARAA[i] = tmp$varAA
}
df = data.frame(Cycle=rep(0:100,4),
                Source=rep(c("Vg","Va","Vd","Vaa"),each=101),
                Variance=c(VARG,VARA,VARD,VARAA))

ggplot(df,aes(x=Cycle,y=Variance,color=Source))+
  geom_line(size=1)+
  guides(alpha=FALSE)+
  theme_bw()

df2 = data.frame(Cycle=0:100,
                 Mean=MEAN)
ggplot(df2,aes(x=Cycle,y=Mean))+
  geom_line(size=1)+
  guides(alpha=FALSE)+
  theme_bw()
