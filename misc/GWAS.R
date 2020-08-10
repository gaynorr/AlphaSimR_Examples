#Example GWAS analysis for Q+K models
library(AlphaSimR) 
library(rrBLUP) #For K and Q+K GWAS
library(qqman) #For Manhattan and Q-Q plots

#Global parameters
nQtlPerChr = 10
nSnpPerChr = 500 #must be > nQtlPerChr

#Simulate mapping population
FOUNDERPOP = runMacs(nInd=1000,
                     nChr=10,
                     segSites=nSnpPerChr,
                     split=20)

SP = SimParam$
  new(FOUNDERPOP)$
  restrSegSites(minSnpPerChr=nSnpPerChr,
                minQtlPerChr=nSnpPerChr,
                overlap=TRUE)$
  addTraitA(nQtlPerChr,gamma=TRUE)$
  addSnpChip(nSnpPerChr)$
  setVarE(H2=0.8)

pop = newPop(FOUNDERPOP)

#Create 2 sub populations
popA = pop[1:500]
popB = pop[501:1000]
for(i in 1:6){
  popA = selectCross(popA,
                     nInd=10,
                     nCrosses=50,
                     nProgeny=10)
  if(i<4){
    popB = selectCross(popB,
                       nInd=10,
                       nCrosses=50,
                       nProgeny=10)
  }
}
pop = c(popA,popB)


#Format genotypes for rrBLUP
rawGeno = pullSnpGeno(pop)
freq = colMeans(rawGeno)/2
MAF = apply(cbind(freq,1-freq),1,min)
rawGeno = rawGeno-1
geno = as.data.frame(t(rawGeno))
geno = data.frame(snp=row.names(geno),
                  chr=rep(1:10,each=nSnpPerChr),
                  pos=rep(1:nSnpPerChr,10),
                  geno)
#Create "pheno" data.frames for rrBLUP
pheno = data.frame(gid=names(geno)[-(1:3)],
                   trait1=pop@pheno[,1])
phenoQ = data.frame(gid=pheno$gid,
                    subPop=factor(rep(c("a","b"),each=500)),
                    trait1=pheno$trait1)
#Find QTL locations within SNP chip
qtl = paste0("SNP_",
             SP$traits[[1]]@lociLoc+
               rep(0:9*nSnpPerChr,each=nQtlPerChr))
##GWAS
#No adjustments----
model0 = data.frame(snp=geno$snp,
                    chr=geno$chr,
                    pos=geno$pos,
                    trait1=rep(NA_real_,10*nSnpPerChr))
for(i in 1:(10*nSnpPerChr)){
  if(MAF[i]>=0.05){
    #Fit linear model with aov and extract p-value
    tmp = summary(aov(pheno$trait1~rawGeno[,i]))[[1]][1,5]
    #Check for bad p-values
    #Check for markers confounded with structure
    if(is.na(tmp)){ 
      model0$trait1[i] = 1
    }else{
      model0$trait1[i] = tmp
    }
  }else{
    model0$trait1[i] = 1
  }
}
#Account for p-value=0
model0$trait1[model0$trait1==0] = 1e-300

#Adjust for population structure (Q)----
modelQ = model0
for(i in 1:(10*nSnpPerChr)){
  if(MAF[i]>=0.05){
    #Fit linear model with aov and extract p-value
    tmp = summary(aov(pheno$trait1~phenoQ$subPop+rawGeno[,i]))
    tmp = tmp[[1]][2,5]
    #Check for markers confounded with structure
    if(is.na(tmp)){ 
      modelQ$trait1[i] = 1
    }else{
      modelQ$trait1[i] = tmp
    }
  }else{
    modelQ$trait1[i] = 1
  }
}
#Account for p-value=0
modelQ$trait1[modelQ$trait1==0] = 1e-300

#Adjust for kinship (K)----
modelK = GWAS(pheno,geno,plot=FALSE)
modelK$trait1 = 10^(-modelK$trait1) #Revert to p-value

#Adjust for structure and kinship (Q+K)----
modelQK = GWAS(phenoQ,geno,fixed="subPop",plot=FALSE)
modelQK$trait1 = 10^(-modelQK$trait1) #Revert to p-value

#Manhattan plots
tiff("manhattan.tiff",height=480,width=480*1.5)
op = par(mfrow=c(2,2),mai=c(.9,.9,0.9,0.9),
         mar=c(2.5,2.5,1,1)+0.1,
         mgp = c(1.5,0.5,0))
manhattan(model0,chr="chr",bp="pos",p="trait1",snp="snp",
          highlight=qtl,main="Unadjusted")
manhattan(modelQ,chr="chr",bp="pos",p="trait1",snp="snp",
          highlight=qtl,main="Adjusted for Q")
manhattan(modelK,chr="chr",bp="pos",p="trait1",snp="snp",
          highlight=qtl,main="Adjusted for K")
manhattan(modelQK,chr="chr",bp="pos",p="trait1",snp="snp",
          highlight=qtl,main="Adjusted for Q+K")
par(op)
dev.off()
system(paste("open","manhattan.tiff"))

#Q-Q plots
tiff("qq.tiff",height=480,width=480*1.5)
op = par(mfrow=c(2,2),mai=c(.9,.9,0.9,0.9),
         mar=c(2.5,2.5,1,1)+0.1,
         mgp = c(1.5,0.5,0))
qq(model0$trait1,main="Unadjusted")
qq(modelQ$trait1,main="Adjusted for Q")
qq(modelK$trait1,main="Adjusted for K")
qq(modelQK$trait1,main="Adjusted for Q+K")
par(op)
dev.off()
system(paste("open","qq.tiff"))

