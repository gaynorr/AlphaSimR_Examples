# Using msprime (a Python library) to replace MaCS

library(AlphaSimR)
library(reticulate) # Requires Python

# py_install("msprime") # Install msprime library

msprime = import("msprime")

# Run msprime simulation
tree_sequence = msprime$
  simulate(
    sample_size=20L, 
    Ne=100L, 
    length=1e8, 
    recombination_rate=1e-8, 
    mutation_rate=2.5e-8)

# Extract genetic map position
sites = tree_sequence$tables$sites$asdict()
pos = sites$position * 1e-8 # Convert to Morgans
pos = pos - pos[1] # Set first position to zero

# Extract haplotypes
haplo = t(tree_sequence$genotype_matrix())

# Create an AlphaSimR founder population
founderPop = newMapPop(genMap=list(pos),
                       haplotypes=list(haplo))
