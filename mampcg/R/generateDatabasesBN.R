
#'##############################################################
#' function for generating a random database from net
#' arguments:
#' @param net: net to sample from (object of bn class (bnlearn))
#' @param numberFiles: number of variants to generate
#' @param numberSamples: number of samples to generate
#' @param path: base path where databases must be stored
#' @param filename: path to the net file. This path is used to extract
#'             the name of the net and to compose the file name 
#'             where to store the database
#'##############################################################
generateDatabases <- function(net, numberFiles, numberSamples, path, filename){
  # gets the file name without extension
  base.filename <- strsplit(filename, '[.]')[[1]][1]
  
  # generates teh required number of databases
  for(i in 1:numberFiles){
    cat(" ... generation for variant: ",i,"....\n")
    # generate the database sampling from 
    db <- cpdist(net, nodes=nodes(net), evidence=TRUE, cluster = NULL, 
                 method = "ls", n=numberSamples, debug = FALSE)
    
    # stores in the corresponding file
    dbname <- paste(base.filename,numberSamples,sep="-")
    dbname <- paste(dbname,i,sep="-")
    dbname <- paste(dbname,".db",sep="")
    dbname <- paste(path,dbname,sep="")
    write.table(db,file=dbname,col.names=TRUE, row.names=FALSE,sep=",")
  }
}

#'##############################################################
#' method receving the name name of the net and generating
#' the corresponding set of databases
#' arrguments:
#' @param netFileName number of the net to generate samples from
#' @param numberFiles number of variants to generate
#' @param numberSamples number of samples to store in each variant
#' @param pathnet path where net files are stored
#' @param pathdb path where databases must be stored
#'##############################################################
generateDatabasesForNet <- function(netFileName, numberFiles=1, numberSamples, 
                                    pathnet, pathdb){
  # composes the db name removing the extension ".net"
  # gets the file name without extension
  cat("Net: ",netFileName," samples: ",numberSamples," numberFiles: ", numberFiles, "\n")
  baseName <- strsplit(netFileName, '[.]')[[1]][1]
  dbFileName<- paste(baseName,".db",sep="")
  
  # composes path names for net and db
  netPathName <- paste(pathnet,netFileName,sep="")
  dbPathname <- paste(pathdb,baseName,sep="")
  dbPathname <- paste(dbPathname,"/",sep="")
  if (!file.exists(dbPathname)){
    dir.create(dbPathname)
  }
  
  # into it creates a new subfolder for the number of samples
  dbPathname <- paste(dbPathname,numberSamples,sep="")
  dbPathname <- paste(dbPathname,"/",sep="")
  if (!file.exists(dbPathname)){
    dir.create(dbPathname)
  }

  # reads the net: this method is defined in dataBaseGeneration.R file
  net <- readNetFile(netPathName,debug=FALSE)
  
  # generates the databases: this method is defined in dataNaseGeneration.R file
  generateDatabases(net,numberFiles,numberSamples,dbPathname,dbFileName)
}

#'##############################################################
#' generate data bases for all the nets stored in the given path
#' arguments:
#' @param pathnet path to folder where nets are stored
#' @param pathdb path where databases must be stored
#' @param numberSamples number of samples to generate
#' @param numberFiles number of databases to generate
#'##############################################################
generateDatabasesForNets <- function(pathnet, pathdb, numberSamples, numberFiles){
  # gets the complete list of net files stored in pathnet folder
  nets <- list.files(pathnet,pattern="*.net",include.dirs=FALSE)

  # generates data bases for every net
  sapply(nets,generateDataBasesForNet, numberFiles=numberFiles,
         numberSamples, pathnet=pathnet,pathdb=pathdb)
}

#'##############################################################
#' generate data bases for a given net (with netname)
#' arguments:
#' @param pathnet path where nets are stored
#' @param netname name of the Bayesian network to sample from
#' @param pathdb path where databases must be stored
#' @param numberSamples number of samples to generate
#' @param numberFiles number of databaes to generate
#'##############################################################
generateDatabasesForNetName <- function(pathnet,netname,pathdb,numberSamples, numberFiles){
  nets <- list.files(pathnet,pattern=paste0(netname,"*.net"),include.dirs=FALSE)

  # generates data bases for every net
  sapply(nets,generateDatabasesForNet, numberFiles,
         numberSamples=numberSamples, pathnet=pathnet,pathdb=pathdb)
}


