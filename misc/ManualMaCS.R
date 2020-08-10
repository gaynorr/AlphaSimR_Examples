# This script shows how supply a manual MaCS command to simulate 10 cattle breeds

library(AlphaSimR)

macsCommand = function(size){
  paste(1.00E+08,"-t",9.00E-06,"-r",3.60E-06,"-I 10",paste(rep(size/10,10),collapse=" "),
        "-eN 0.011 1.33 -eN 0.019 2.78 -eN 0.036 3.89 -eN 0.053 11.11 -eN 0.069 16.67 -eN 0.431 22.22 -eN 1.264 27.78 -eN 1.819 38.89 -eN 4.875 77.78 -eN 6.542 111.11 -eN 9.319 188.89 -eN 92.097 688.89 -ej 0.136 2 1 -ej 0.137 3 1 -ej 0.138 5 4 -ej 0.139 6 4 -ej 0.140 8 7 -ej 0.141 9 7 -ej 0.142 10 7 -ej 1.111 4 1 -ej 1.112 7 1 -s",sample.int(1e8,1))
}

nInd = 1000 #Each breed will be 1/10th this size (must be divisible by 10)
nChr = 10
nSnpPerChr = 200
nQtlPerChr = 100

founderPop = runMacs(nInd,nChr,nQtlPerChr+nSnpPerChr,
                     manualCommand=macsCommand(2*nInd), #Command for simulating haplotypes
                     manualGenLen=1) #Genetic length of chromosomes

SP = SimParam$new(founderPop[1:(nInd/10)]) #Use the first breed as the founder population

#Save each breed as an element in a list
breeds = vector("list",10)
for(i in 1:10){
  breeds[[i]] = newPop(founderPop[((i-1)*(nInd/10)+1):(i*nInd/10)])
}

