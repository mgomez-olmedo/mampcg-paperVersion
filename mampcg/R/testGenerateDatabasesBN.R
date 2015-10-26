
# sets relevant paths to nets folder and databases
path <- "./networks/"
pathddbb <- "./ddbb/"

# generates the databases for all these sample sizes
ns <- c(500,1000,5000,10000,50000)

# sets the number of databases to generate for each dataset
numberFiles <- 30

# sets the name of the net used for generating samples
netname <- "alarm"

# considers each sample size
for(i in 1:length(ns)){
  # selects the corresponding sample size
  nsamples <- ns[i]
  
  # generate the databases
  generateDatabasesForNetName(path, netname, pathddbb, nsamples, numberFiles)
}
