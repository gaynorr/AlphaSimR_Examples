# This is the first script of the simulation
# All other scripts are called by this script
library(AlphaSimR)

# Load global parameters
source("GlobalParameters.R")

nReps = 10 #Number of replications for experiment
burninYears = 20
futureYears = 20
# Initialize variables for results
hybridCorr = inbredMean = hybridMean = inbredVar = hybridVar =
  rep(NA_real_,burninYears+futureYears)
output = list(inbredMean=NULL,
              inbredVar=NULL,
              hybridMean=NULL,
              hybridVar=NULL,
              hybridCorr=NULL)
# Save results
saveRDS(output,"Scenario1.rds")
saveRDS(output,"Scenario2.rds")

# Loop through replications
for(rep in 1:nReps){
  # Create initial parents and set testers and hybrid parents
  source("CreateParents.R")

  # Fill breeding pipeline with unique individuals from initial parents
  source("FillPipeline.R")

  # p-values for GxY effects
  P = runif(burninYears+futureYears)

  # Cycle years
  for(year in 1:burninYears){ #Change to any number of desired years
    cat("Working on year:",year,"\n")
    p = P[year]
    source("UpdateParents.R") #Pick new parents based on last year's data
    source("UpdateTesters.R") #Pick new testers and hybrid parents
    source("AdvanceYear.R") #Advances yield trials by a year
    source("UpdateResults.R") #Track summary data
  }

  #Save burn-in to load later for scenario 2
  save.image("tmp.rda")

  #Scenario 1
  cat("Working on Scenario 1\n")
  for(year in (burninYears+1):(burninYears+futureYears)){
    cat("Working on year:",year,"\n")
    p = P[year]
    source("UpdateParents.R") #Pick new parents based on last year's data
    source("UpdateTesters.R") #Pick new testers and hybrid parents
    source("AdvanceYear.R") #Advances yield trials by a year
    source("UpdateResults.R") #Track summary data
  }

  # Report results for scenario 1
  output = readRDS("Scenario1.rds")
  output = list(inbredMean=rbind(output$inbredMean,inbredMean),
                inbredVar=rbind(output$inbredVar,inbredVar),
                hybridMean=rbind(output$hybridMean,hybridMean),
                hybridVar=rbind(output$hybridVar,hybridVar),
                hybridCorr=rbind(output$hybridCorr,hybridCorr))
  saveRDS(output,"Scenario1.rds")

  # Scenario 2
  load("tmp.rda")
  cat("Working on Scenario 2\n")
  for(year in (burninYears+1):(burninYears+futureYears)){
    cat("Working on year:",year,"\n")
    p = P[year]
    source("UpdateParents.R") #Pick new parents based on last year's data
    source("UpdateTesters.R") #Pick new testers and hybrid parents
    source("AdvanceYearAlt.R") #Advances yield trials by a year
    source("UpdateResults.R") #Track summary data
  }

  # Report results for scenario 2
  output = readRDS("Scenario2.rds")
  output = list(inbredMean=rbind(output$inbredMean,inbredMean),
                inbredVar=rbind(output$inbredVar,inbredVar),
                hybridMean=rbind(output$hybridMean,hybridMean),
                hybridVar=rbind(output$hybridVar,hybridVar),
                hybridCorr=rbind(output$hybridCorr,hybridCorr))
  saveRDS(output,"Scenario2.rds")
}

#Delete tmp file
file.remove("tmp.rda")
