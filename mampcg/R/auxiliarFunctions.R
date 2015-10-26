
library(bnlearn)

#'##############################################################
#' function for reading the data of the net from a file
#' arguments:
#' @param fileName complete path to the file to read. The function
#'             read.net is defined in bnlearn library
#'##############################################################
readNetFile <- function(fileName,debug=FALSE){
  cat("Filename: ",fileName, "\n")
  net <- bnlearn::read.net(fileName,debug)
}

#'##############################################################
#' function for reading a database file
#' arguments:
#' @param filedb path to the database to read
#' @param netExt flag to show if data must be considered as
#'        factor or not
#'##############################################################
readDataBase <- function(filedb, netExt){
  db <- read.csv(filedb,header=TRUE)
  colnames <- colnames(db)
  
  # converts everything into factor if netExt is "net": then
  # data come from discrete Bayesian networks
  if (netExt == "net"){
    db[,colnames] <- lapply(db[,colnames] , factor)
  }
  
  # return db
  return(db)
}
