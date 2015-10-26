# definition of global parameters: paths, base file name (artificial),
# sample sizes, vector of edges probs for the variants to consider
# (1. 20% for directed edges, 20% undirected, 60% bidirected;
#  2. 20% directed, 60% undirected, 20% bidirected;
#  3. 60% directed, 20% undirected, 20% bidirected)
# number of nets to generate, number of nodes (15), number of databases
# to generate during the sample procedure and seed (NOTE: remove this
# sentence for a real random behaviour; in this cases is fixed just to
# focus the software on the tree models used for the experiments described
# in the paper)
pathDb <- "./ddbb/"
pathNet <- "./networks/"
baseFileName <- "artificial"
#sampleSizes <- c(500, 1000, 5000, 10000, 50000)
sampleSizes <- c(200)

# these are the probs for directed, undirected and bidirected used for
# artificial1, artificial2 and artificial3 respectively
edgesProbs <- list(c(0.2, 0.2, 0.6), c(0.2, 0.6, 0.2), c(0.6, 0.2, 0.2))
numberNets <- 3
numberNodes <- 15
variants <- 5

# This parameter defines the threshold used for the evidence. We have
# used 0.2 for artificial1 and artificial3 and 0.9 for artificial2
threshold <- 0.2

# use parallel execution if possible. The argument defines the number of
# cores to use
cl=makeCluster(2)

# normal procedure
# 1. generates a new graphical model
graphModel <- generateRandomMAMPCG(numberNodes, edgesProbs[[1]])

# 2. prepare the model for sampling. This also requires to get a valid
# name for the model (it will be artificial, bit the id will be assigned
# examining the folder considered for storing the models (pathNet))
completeModel <- prepareMAMPCGForSampling(graphModel, baseFileName, pathNet)

# 3. If needed, store the model for posterior use. The procedure
# will check all the files named artificial and will add a create
# a new one. This method only saves complete models defined by
# edges, a bnet and a set of inmoralities
storeModel(completeModel,pathNet)

# 4. The model can be used to generate the samples
generateDatabases(completeModel, variants, sampleSizes, threshold, pathDb, cl)

# 5. We can retrieve a previously generated model for producing new datasets
completeModel <- retrieveModel("artificial1",pathNet)

# 6. Now it is possible to generateDatasets as before....
# generateDatabases(completeModel, variants, sampleSizes, threshold, pathDb, cl)