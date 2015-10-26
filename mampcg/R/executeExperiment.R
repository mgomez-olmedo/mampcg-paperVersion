
#'##############################################################
#' learning function. It prepares the creation of the object of 
#' MampcSearch and once create calls the method for running the
#' learning algorithm
#' arguments
#' @param pathnet path to nets
#' @param pathdb path to ddbb
#' @param netName net name
#' @param dbId ddbb id to learn from
#' @param numberSamples number of samples of the ddbb to analyze
#' @param moral way to consider v-structures
#' @param pc flag to show the algorithm to use (mampc or pc)
#' @param netExt specify the origin of the true net (bn or mamp graph)
#' @param debug flag to activate/deactivate the use of debug messages
#'##############################################################
learn <- function(pathnet, pathdb, netName, dbId, numberSamples, moral, pc, netExt, debug){
  # sets edges and net to null
  edges <- NULL
  net <- NULL
  
  # the true net is a a real BN
  if (netExt == "net"){
    # compose the name of the net
    netFileName <- paste(pathnet,netName,sep="")
    netFileName <- paste(netFileName,".net",sep="")
    
    # reads the net
    net <- read.net(netFileName,FALSE)
    if (debug){
      cat("Read information of true net..................\n")
    }
  }
  else{
    # read edges from a RDS file (it is a mamp graph and can not be 
    # described with a BN)
    edgesFileName <- paste(pathnet, netName,sep="")
    edgesFileName <- paste(edgesFileName, ".mampcg", sep="")
    
    # read the RDS file with model information
    model <- readRDS(edgesFileName)
    if (debug){
      cat("Read information about edges: number of nodes = ",
                                          length(bnlearn::nodes(model$bnet)),"\n")
      cat("number of edges = ",nrow(model$edges),"\n")
    }
    
    # sets net and edges
    net <- model$bnet 
    edges <- model$edges
  }
  
  # compose the name of the database to analyze
  # As an example: alarm-1000-10.db it would be the name
  # of a database generated from alarm netwotk, with 1000 
  # samples and being the variant number 10 for this net 
  # and number of samples
  dbName <- paste(netName,"-",sep="")
  dbName <- paste(dbName,numberSamples,sep="")
  dbName <- paste(dbName,"-",sep="")
  dbName <- paste(dbName,dbId,sep="")
  dbName <- paste(dbName,".db",sep="")
  dbCompleteName <- paste(pathdb,dbName,sep="")
  
  # reads the data base to analyze
  db <- readDataBase(dbCompleteName, netExt)
  if (debug){
    cat("Read ddbb: ",dbCompleteName,"\n")
  }
  
  # now it is time to learn. For this purpose is required to create
  # an object on MampcSearch, passign as arguments the net, data base,
  # moral flag, pc flag, set of edges (null if the initial net is a BN)
  # and debug flag
  bInfo <- buildObject(net, db, moral=moral, pc=pc, edges=edges, debug=debug)
  
  # calls the method devoted to execute the learning algorithm
  bInfo$learn()
  
  # return bInfo
  return(bInfo)
}

#'##############################################################
#' This function defines the execution of the learning algorithm
#' arguments:
#' @param pathnet: path to nets
#' @param pathdbsample: path to ddbbs
#' @param netName: net name
#' @param dbId: variant to consider (each ddbb is identified with a different id)
#' @param samples: number of samples of the ddbb to consider
#' @param mode: mode of composing the true net (from a BN defined in a net file
#'        of from a rds file with the definition of the edges in case of a mamp
#'        graph) 
#' @param debug: flag to show if trace messages must be used or not
#'##############################################################
execute <- function(pathnet, pathdbsample, netName, dbId, samples, mode, debug){
  # sets model and modelPC (the results of the learning algorithm
  # to NULL)
  model <- NULL
  modelPC <- NULL

  if (debug){
    cat("Learning from variant: ",dbId,"\n")
  }

  # learn the model with mampc algorithm. By default debug mode is. The basic
  # function learn is defined in the file learn.R
  cat("Sample size: ",samples, " - variant: ", dbId, "\n");
  cat("MAMPCG result\n");
  try(model <- learn(pathnet,pathdbsample,netName,dbId,samples,moral=FALSE, 
                     pc=FALSE, mode, debug))
  
  if (debug){
    cat ("   learnt with mamp algorithm\n")
  }

  # learn the model with pc algorithm: only if mode is net (learning from 
  # Bayesian networks)
  if (mode == "net"){
    cat("PC result\n");
    try(modelPC <- learn(pathnet,pathdbsample,netName,dbId,samples,moral=FALSE, 
                         pc=TRUE,mode,debug))
  
    if (debug){
      cat("    learnt with pc algorithm\n")
    }
  }

  # compare learnt models against the true model
  # first compare the model learnt with mampc algorithm
  if (!is.null(model)){
    vs <- model$compareStructure()
  }
  else{
    # if the algorithm was no able to finish, considers NA as values
    # for recall, precision, recallVs and precisionVs
    vs <- list(recall=NA, precision=NA, recallVs=NA, precisionVs=NA)
  }

  # store the results of the comparison: this method is defined in utilResults.R
  # script
  storeNetResults(netName, "mampcg", samples, dbId, c(vs$recall,vs$precision,
                                                      vs$recallVs, vs$precisionVs))

  # now compares result pf pc algorithm and true net
  if (!is.null(modelPC)){
    vsPC <- modelPC$compareStructure()
  }
  else{
    # if the algorithm was no able to finish, considers NA as values
    # for recall, precision, recallVs and precisionVs
    vsPC <- list(recall=NA, precision=NA, recallVs=NA, precisionVs=NA)
  }

  # store the results of the comparison: this method is defined in utilResults.R
  # script
  storeNetResults(netName, "pc", samples, dbId, c(vsPC$recall,vsPC$precision,
                                                  vsPC$recallVs, vsPC$precisionVs))

  # return the results of the comparison
  return(c(vs$recall, vs$precision, vs$recallVs, vs$precisionVs, 
           vsPC$recall, vsPC$precision, vsPC$recallVs, vsPC$precisionVs))
}
