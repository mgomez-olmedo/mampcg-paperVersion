
#'##############################################################
#' stores the data of a learning process on netName with
#' a certain method, number of samples, id of data base file
#' arguments:
#' @param netName name of the net. The results will be related to the
#'            folder belonging to this net
#' @param method this flag shows the algorithm producing the results
#'           (mamp - pc)
#' @param samples number of samples of the database used for learning
#'           (this name is used to organize the results into folders)
#' @param id id of the database. The is will be included in the filename
#'            including the results
#' @param results results to store
#'##############################################################
storeNetResults <- function(netName, method, samples, id, results){
  # checks if the folder where to store the results exist
  basicPath <- "./results/"
  if (!file.exists(basicPath)){
    dir.create(basicPath)
  }
  
  # checks if the folder for the net exist
  path <- paste(basicPath,netName,sep="/")
  if (!file.exists(path)){
    dir.create(path)    
  }
  basicPath <- path
  
  # checks if the folder for the method exist
  path <- paste(basicPath,method,sep="/")
  if (!file.exists(path)){
    dir.create(path)    
  }
  basicPath <- path
  
  # checks the existence of the folder with the number of samples
  path <- paste(basicPath,samples,sep="/")
  if (!file.exists(path)){
    dir.create(path)    
  }
  
  # write the file
  fileName <- paste(netName,"-",sep="")
  fileName <- paste(fileName,samples,sep="")
  fileName <- paste(fileName,"-",sep="")
  fileName <- paste(fileName,id,sep="")
  fileName <- paste(fileName,".dat",sep="")
  finalPath <- paste(path,fileName,sep="/")
  write.table(results, file = finalPath, row.names=FALSE, col.names=FALSE, sep=",")
}

#'##############################################################
#' method to load the results of a certain execution
#' arguments: 
#' @param netName name of the net of interest
#' @param method method which produced the results
#' @param samples number of samples of the database
#' @param id id of the database
#' NOTE: this method geneates the results in a similar way as
#' they appear in the paper
#'##############################################################
loadNetResults <- function(netName, method, samples, id){
  # compone el nombre del path
  basicPath <- "./results"
  path <- paste(basicPath,netName,sep="/")
  path <- paste(path,method,sep="/")
  path <- paste(path,samples,sep="/")
  path <- paste(path,netName,sep="/")
  path <- paste(path,samples,sep="-")
  path <- paste(path,id,sep="-")
  path <- paste(path,".dat",sep="")

  # se lee el archivo
  data <- as.data.frame(c(0,0,0,0))
  if (file.exists(path)){
     data <- read.csv(path,header=FALSE)
  }
  return(data)
}

#'##############################################################
#' method to load the results of a certain execution
#' arguments: 
#' @param netName name of the net of interest
#' @param method method which produced the results
#' @param samples number of samples of the database
#' @param repetitions variants for each sample size
#'##############################################################
processResults <- function(netName, method, samples, repetitions){
  
  # se crea la matriz donde se almacenaran al final los resultados
  resultsMean <- matrix(0,4,length(samples))
  resultsDev <- matrix(0,4,length(samples))
  
  # considera cada posible tamaño de conjunto de datos
  for(i in 1:length(samples)){
    size <- samples[i]
    
    # se crea matriz de resultados para cada tamaño muestral
    partialResults <- matrix(0,4,repetitions)
    
    for(j in 1:repetitions){
      # se leen los datos de cada conjunto de datos
      result <- loadNetResults(netName, method, size, j)
      
      # se vuelcan sobre la matriz
      partialResults[,j] <- result[,1]
    }
    
    cat("Partial results............................\n")
    print(partialResults)
    cat("-------------------------------------------\n")
    
    # now computes mean and deviation for this results
    for(j in 1:4){
      vals <- partialResults[j,]
      
      # counts all meaningful values (non cero, non NA
      vals <- vals[!is.na(vals)]
      
      # compute mean and deviation
      resultsMean[j,i] <- mean(vals)
      resultsDev[j,i] <- sd(vals)
    }
  }
  
  # return results
  return(list(mean=resultsMean,dev=resultsDev))
}


#'##############################################################
#' generate a latex table reading the data pased as third argument
#' arguments: 
#' @param netName name of the net
#' @param samples vector of sample sizes
#' @param results results of the execution to put into the table
#'##############################################################
generateLatexTableFromData <- function(netName, samples, results){
  
  # se crea la matriz donde se almacenaran al final los resultados
  resultsMean <- matrix(0,8,length(samples))
  resultsDev <- matrix(0,8,length(samples))

  for(i in 1:length(samples)){

    # compute mean and deviations removing NAs
    resultsMean[,i] <- apply(results[[i]],1,function(x){
      x <- x[!is.na(x)]
      mean(x)
    })
  
    resultsDev[,i] <- apply(results[[i]],1,function(x){
      x <- x[!is.na(x)]
      sd(x)
    })
  }
  
  sink()
  sink()
  options(digits=2)
  traza <- paste(netName,"-results",sep="")
  sink(traza)
  options(digits=2)
  
  # genera la tabla latex con los resultados
  cat("\\begin{table}[h!]\n")
  cat("\\centering\n")
  cat("\\begin{tabular}{|c|c|c|c|c|c|c|}\\hline")
  cat(" size & &")
  for(i in 1:length(samples)){
    cat(samples[i], " ")
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\hline\\hline\n")
  
  cat("\\multirow{4}{*}{Our algorithm}\n")
  
  # recall data
  cat("  & RA & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    cat(resultsMean[1,i], " $\\pm$ ", resultsDev[1,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{2-7}\n")
  
  # precision data
  cat("  & PA & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    cat(resultsMean[2,i], " $\\pm$ ", resultsDev[2,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{2-7}\n")
  
  # recallVS data
  cat("  & RT & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    if (is.na(resultsMean[3,i])){
      resultsMean[3,i] <- 0
      resultsDev[3,i] <- 0
    }
    cat(resultsMean[3,i], " $\\pm$ ", resultsDev[3,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{2-7}\n")
  
  # prevVS data
  cat("  & PT & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    if (is.na(resultsMean[4,i])){
      resultsMean[4,i] <- 0
      resultsDev[4,i] <- 0
    }
    cat(resultsMean[4,i], " $\\pm$ ", resultsDev[4,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{1-7}\n\\hline\\hline\n")
  
  cat("\\multirow{4}{*}{Meek's algorithm}\n")
  
  # recall data
  cat("  & RA & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    cat(resultsMean[5,i], " $\\pm$ ", resultsDev[5,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{2-7}\n")
  
  # precision data
  cat("  & PA & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    cat(resultsMean[6,i], " $\\pm$ ", resultsDev[6,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{2-7}\n")
  
  # recallVS data
  cat("  & RT & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    if (is.na(resultsMean[7,i])){
      resultsMean[7,i] <- 0
      resultsDev[7,i] <- 0
    }
    cat(resultsMean[7,i], " $\\pm$ ", resultsDev[7,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{2-7}\n")
  
  # prevVS data
  cat("  & PT & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    if (is.na(resultsMean[8,i])){
      resultsMean[8,i] <- 0
      resultsDev[8,i] <- 0
    }
    cat(resultsMean[8,i], " $\\pm$ ", resultsDev[8,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{1-7}\n")
  cat("\\end{tabular}\n")
  cat("\\end{table}\n")
  sink()
  sink()
  sink()
}

#'##############################################################
#' generate a latex table reading the data from the result files
#' arguments: 
#' @param netName name of the net
#' @param samples vector of sample sizes
#' @param repetitions number of databases (variants) to consider
#'##############################################################
generateLatexTableFromFiles <- function(netName, samples, repetitions){
  # gets the data from the files and process them. This is done by
  # the function processResults defined in processResults.R file
  resultsMampcg <- processResults(netName, "mampcg", samples, repetitions)
  resultsPc <- processResults(netName, "pc", samples, repetitions)

  # generate the table
  sink()
  sink()
  options(digits=2)
  traza <- paste(netName,"-results2",sep="")
  sink(traza)
  options(digits=2)
  
  # genera la tabla latex con los resultados
  cat("\\begin{table}[h!]\n")
  cat("\\centering\n")
  cat("\\begin{tabular}{|c|c|c|c|c|c|c|}\\hline")
  cat(" size &  &")
  for(i in 1:length(samples)){
    cat(samples[i], " ")
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\hline\\hline\n")
  
  cat("\\multirow{4}{*}{Our algorithm}\n")
  
  # recall data
  cat("  & RA & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    cat(resultsMampcg$mean[1,i], " $\\pm$ ", resultsMampcg$dev[1,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{2-7}\n")
  
  # precision data
  cat("  & PA & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    cat(resultsMampcg$mean[2,i], " $\\pm$ ", resultsMampcg$dev[2,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{2-7}\n")
  
  # recallVS data
  cat("  & RT & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    cat(resultsMampcg$mean[3,i], " $\\pm$ ", resultsMampcg$dev[3,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{2-7}\n")
  
  # prevVS data
  cat("  & PT & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    cat(resultsMampcg$mean[4,i], " $\\pm$ ", resultsMampcg$dev[4,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{1-7}\n\\hline\\hline\n")
  
  cat("\\multirow{4}{*}{Meek's algorithm}\n")
  
  # recall data
  cat("  & RA & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    cat(resultsPc$mean[1,i], " $\\pm$ ", resultsPc$dev[1,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{2-7}\n")
  
  # precision data
  cat("  & PA & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    cat(resultsPc$mean[2,i], " $\\pm$ ", resultsPc$dev[2,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{2-7}\n")
  
  # recallVS data
  cat("  & RT & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    cat(resultsPc$mean[3,i], " $\\pm$ ", resultsPc$dev[3,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{2-7}\n")
  
  # prevVS data
  cat("  & PT & ")
  
  # salida de valores de recall para cada tamaño muestral
  for(i in 1:length(samples)){
    cat(resultsPc$mean[4,i], " $\\pm$ ", resultsPc$dev[4,i])
    if (i != length(samples)){
      cat(" & ")
    }
  }
  cat("\\\\\\cline{1-7}\n")
  cat("\\end{tabular}\n")
  cat("\\end{table}\n")
  sink()
  sink()
  sink()
}

# this is an example about the way of retrieving the result of
# a given execution
#data <- loadNetResults("alarm","pc",100,1)

