source("R/significanceTest.R")

# these sentences show how to perform the analysis for a given
# network. The experiments for this network must be already
# available
netName <- "asia"
samples <- c(500, 1000, 5000, 10000, 50000)
repetitions <- 30

# prepares the output to a file
sink()
sink()
sink()

# compose the name of the file to generate with the results
traza <- paste0(netName,"-tests")

# redirect output to the file
sink(traza)

# perform the test
results <- significanceTest(netName, samples, repetitions)
sink()
sink()
