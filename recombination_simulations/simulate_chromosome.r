## Simulate a pig population and calculate number of recombinations. 

# Load required packaged. hash, readr and AlphaSimR are availible on CRAN
library(AlphaSimR)
library(readr)
library(hash)


options(scipen=999)
args = commandArgs(trailingOnly=TRUE)

## Command line options:
## This script works on a single chromosome at a time. The chromosome number is given as a command line option.
chr = as.numeric(args[1])

## Input files. We take in a real genotype file and a pedigree file.
## The pedigree file is used to simulate the full pedigree (using new ancestral haplotypes).
## The genotype file is used to determine the number of segregating sites to simulate, and to mask the output genotypes (to take into account multiple chip densities and spontaneous missing markers). 
## read_table2 is taken from the readr file and allows fairly large files to be read-in.

real_genotypes <- data.frame(read_table2(paste0("genotypes.", chr, ".txt"), col_names = FALSE))
pedigree_initial <- read.table("pedigree.txt", header = FALSE, stringsAsFactors = FALSE, colClasses = "character")


## Simulation parameters
## Number of founders for AlphaSimR to simulate. This value just depends on the pedigree and not the genotype file (or chromosome).
## Can set this manually if AlphaSimR throws an error message.  

n_founders = sum(pedigree_initial[,2] == 0 | pedigree_initial[,3] == 0)


nSnp = ncol(real_genotypes) - 1 # The first column of the genotype file are ids.

## Run MACS:
## Cattle is a reasonable stand-in for the ancestral history of a livestock species. 
## It's possible to use a custom Macs command as well (or I think to use phased haplotypes for the founders).
print(paste("Data read in, starting coalescent for Chromosome ", chr, ", n_founders =", n_founders))

founderPop <- runMacs(nInd = n_founders,
                      nChr = 1,
                      segSites = nSnp, species = "CATTLE") 
SP <- SimParam$new(founderPop)
SP$setTrackPed(TRUE)
SP$setTrackRec(TRUE)
SP$setGender("no")



## Create a realistic genetic map.
## NOTE: if you don't care about the exact genetic map, most of this section can be skipped and default map length of 100cM will be used.

# These values are estimate map lenghts for pigs taken from Tortereau et al. 2012 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3499283/)
# We simulated sex-specific recombination rates since there are known sex differences in recombination rate in pigs.
table_mat = c(138 ,133 ,130 ,137 ,128 ,170 ,146 ,120 ,144 ,126 ,107 ,117 ,115 ,135 ,115 ,96 ,100 ,87)
table_pat = c(151, 111, 111, 112, 100, 125, 118, 105, 110, 93, 63, 81, 112, 112, 102, 70, 56, 50)

pat_rec = table_pat[chr]
mat_rec = table_mat[chr]

tmp = lapply(SP$genMap,function(x) (pat_rec+mat_rec)/(2*100)*x) # Modifies the recombination rate based on the paternal chromosome.
SP$switchGenMap(tmp)

rec_ratio = mat_rec/pat_rec
SP$setMeiosis(ratio=rec_ratio) # Set's the ratio between maternal and paternal recombinations.


SP$addSnpChip(nSnpPerChr = nSnp)


## Run simulation using the input pedigree.

print(paste("Simulating Chromosome", chr))

pop <- newPop(founderPop)
output <- pedigreeCross(pop,
                        id = pedigree_initial[,1],
                        mother = pedigree_initial[,3],
                        father = pedigree_initial[,2])



## Get the resulting genotypes and pedigree

geno <- pullSnpGeno(output)

pedigree <- data.frame(output@id,
                       output@father,
                       output@mother,
                       stringsAsFactors = FALSE)

masked_genotypes <- geno

## Add a 2% genotyping error. Can remove if errors do not need to be added.

addErrors = function(values){
    newValues = (values + round(runif(length(values), min = 1, max = 2))) %% 3
    return(newValues)

}
error_rate = 0.02
markers = which(rbinom(length(masked_genotypes), 1, error_rate) == 1)
if(length(markers) >0){
    masked_genotypes[markers] = addErrors(masked_genotypes[markers])
}


## Mask genotypes to mimic the input genotype file.

id_to_genotype = hash()
for(i in 1:nrow(pedigree_initial)){
    id_to_genotype[pedigree_initial[i,1]] = -1
}
for(i in 1:nrow(real_genotypes)){
    id_to_genotype[real_genotypes[i,1]] = i
}

masking = as.matrix(real_genotypes[,-1] == 9)

for (i in 1:nrow(masked_genotypes)) {
    
    if (i%%1000 == 0) print(i)

    genotype_id = id_to_genotype[[pedigree_initial[i,1]]]

    # Masks the individual's genotypes to low density.
    # If not in the genotype file, mask all of the genotypes.
    # Otherwise mask using the appropriate genotype.

    if (genotype_id == -1) {
        masked_genotypes[i,] = 9
    }
    else{
        masked_genotypes[i, masking[genotype_id,]] = 9
    }
}



## Calculates and stores the true recombination history

recHist <- SP$recHist
tmp = lapply(1:length(recHist), function(index) {
    ind = recHist[[index]]
    if(!is.null(ind)) {
        pat = ind[[1]][[1]]
        values_pat = nrow(pat) - 1 
        mat = ind[[1]][[2]]
        values_mat = nrow(mat) - 1 
        values_total = values_pat + values_mat
        return(c(index, values_pat, values_mat, values_total))
    }
    return(c(index, NA, NA, NA))
})

recombinations = do.call("rbind", tmp)

recombinations_subset = as.data.frame(recombinations[recombinations[,1] %in% pedigree[,1],])

ped_id_to_real = hash( as.character(pedigree[,1]), pedigree_initial[,1])

rec_ids = sapply(as.character(recombinations_subset[,1]), function(id){
    return(ped_id_to_real[[id]])
})
recombinations_subset[,1] = rec_ids

## Output files.
output_folder = "simulated_data"
dir.create(output_folder)

write.table(pedigree_initial,
            file = file.path(output_folder, "pedigree.txt"),
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

write.table(cbind(pedigree_initial[,1], geno),
            file = file.path(output_folder, paste0("true_genotypes.", chr, ".txt")),
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

write.table(cbind(pedigree_initial[,1], masked_genotypes),
            file = file.path(output_folder, paste0("genotypes.", chr, ".txt")),
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

write.table(recombinations_subset,
            file = file.path(output_folder, paste0("recombinations.", chr, ".txt")),
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)
