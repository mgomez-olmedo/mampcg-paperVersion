
source("R/utilResults.R")

#'##############################################################
#' process results in order to compute the significance tests
#' arguments:
#' @param netName name of the net to analyze
#' @param samples vector with sample sizes
#' @param repetitions number of repetitions for each sample size
#'##############################################################
significanceTest <- function(netName, samples, repetitions){

  # se crea la matriz donde se almacenaran al final los resultados
  resultsTest <- matrix(0,4,length(samples))

  # considera cada posible tamaño de conjunto de datos
  for(i in 1:length(samples)){
    # shows the sample size
    cat("---------------- Sample size: ",samples[i],"----------------------\n")
    size <- samples[i]
  
    # se crea matriz de resultados para cada tamaño muestral
    partialResults1 <- matrix(0,4,repetitions)
    partialResults2 <- matrix(0,4,repetitions)
  
    # gets the data for each variant
    for(j in 1:repetitions){
      # se leen los datos de cada conjunto de datos
      result1 <- loadNetResults(netName, "mampcg", size, j)
      result2 <- loadNetResults(netName, "pc", size, j)
    
      # se vuelcan sobre la matriz
      partialResults1[,j] <- result1[,1]
      partialResults2[,j] <- result2[,1]
    }
    
    # now computes the test for every measure
    for(j in 1:4){
      # shows information about the values to be compared
      cat("Performing tests for measure number ",j,"-----------------------\n")
      
      vals1 <- partialResults1[j,]
      vals2 <- partialResults2[j,]
      
      # shows the serie of values
      cat("Values for mampc: ")
      print(vals1)
      cat("Values for pc: ")
      print(vals2)
      

      # counts all meaningful values (non cero, non NA)
      vals1 <- vals1[!is.na(vals1)]
      vals2 <- vals2[!is.na(vals2)]
      
      cat("Non NA values for mampc: ",length(vals1),"\n")
      cat("Non NA values for pc: ",length(vals2),"\n")
      
      # make the test
      resultTest <- NULL
      try(resultTest <- t.test(vals1, vals2, alt="two.sided",var.equal=TRUE))

      if (!is.null(resultTest)){ 
      	resultsTest[j,i] <- resultTest$p.value
      	
      	# shows test result
      	cat("Performed test --- result: \n")
      	print(resultTest)
      	
      }
      else{
        # shows test could not be performed
        cat("Test can not be computed \n")
      	resultsTest[j,i] <- NA
      }
    }
  }
 
  # return results
  return(resultsTest)
}
