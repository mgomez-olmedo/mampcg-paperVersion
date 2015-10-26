
library(foreach)
library(doParallel)

# This sentences configure the environment for performing the
# experiment
source("R/auxiliarFunctions.R")
source("R/utilResults.R")
source("R/MAMPCGSearch.R")
source("R/executeExperiment.R")

#sets the number of digits for the final report
options(digits=2)

# sets the paths to the databases and to the networks. This has to
# be changed when the software is installed on another machine
pathdb <- "./ddbb/"
pathnet <- "./networks/"

# sets the name of the network to work with
netName <- "barley"

# clean ths sink in order to initiate a new trace if required. If this
# is the case, the last of these three lines must be operative removing
# the comment mark
sink()
sink()
traceFileName <- paste0(netName,"-trace")
sink(traceFileName)

# uses parallelism if possible. The next sentence sets the number of
# cores to use (depending on the execution machine)
registerDoParallel(cores=15)
getDoParWorkers()

# set the different sample sizes to consider
samples <- c(500, 1000, 5000, 10000, 50000)

# sets the number of variants for each sample size
repetitions <- 30

cat("Learning process start\n")

# initializes globalResults with an empty list
globalResults <- list()

# sets this var to show the origin of the true model (net or
# edges). Now it is set  to "net" in order to learn a BN. Change
# this value to edges if the objective if to learn (and compare)
# with respect to a mamp model
trueModel <- "net"

# sets debug flag
debug <- FALSE

# considers each sample size
for(i in 1:length(samples)){
  cat("Learning for sample size: ",samples[i],"\n")
  pathdbsample <- ""
  # composes the path where the ddbb is located
  pathdbsample <- paste(pathdb,netName,sep="")
  pathdbsample <- paste(pathdbsample,samples[i],sep="/")
  pathdbsample <- paste(pathdbsample,"/",sep="")
  cat("ddbb path: ",pathdbsample,"\n")
  
  # creates the result matrix for this sample size: it will
  # contain 8 rows (results of mampc and pc) and as many
  # columns as the number of variants for this sample size
  partialResults <- matrix(NA,8,repetitions)
  
  # considers every repetition with a parallel approach
  partialResults <- foreach(j=1:repetitions, .combine='cbind') %dopar% {
    # learn from the ddbb, with ths corresponding sample size and
    # the variant given by j
    execute(pathnet, pathdbsample, netName, j, samples[i], trueModel, debug)
  }
  
  # stores the results into globalResults
  colnames(partialResults) <- (c(1:repetitions))
  globalResults[[i]] <- partialResults
}

# generate latex table from data (globalResult) (if required). Only
# one of these sentences must be employed. It is included here in
# order to check the algorithm
generateLatexTableFromData(netName, samples, globalResults)

# generate latex table from files (if required)
#generateLatexTableFromFiles(netName, samples, repetitions)

# gets sure no sink is open
sink()
sink()
