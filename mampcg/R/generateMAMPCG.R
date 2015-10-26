library(bnlearn)
library(parallel)


#'###################### SECTION 1: auxiliar functions #########################

#'#############################################################
#' gets a new id for a model analyzing a given folder. If networks
#' folder contains models artificial1....artificial10 the new model
#' will be artificial11
#' arguments:
#' @param baseFileName reference name
#' @param folder to analyze
#' @return new file name
#'#############################################################
getNewFileName <- function(baseFileName, pathNet){
  # gets all the files matching baseFileName in the path
  files <- NULL
  files <- list.files(pathNet, pattern=baseFileName, include.dirs=FALSE)
  
  # gets the biggest id
  maxId <- 0
  if (length(files) != 0){
    for(i in 1:length(files)){
      name <- files[[i]]
      # removes the extension
      basename <- strsplit(name,'[.]')[[1]][1]
      
      # gets the number
      number <- as.numeric(substring(basename,nchar(baseFileName)+1, 
                                     nchar(basename)))
      
      # checks if needed to update maxId
      if (maxId < number){
        maxId <- number
      }
    }
  }
  
  # returns the new name
  name <- paste0(baseFileName,maxId+1)
}

#'#############################################################
#' creates a directed edge between two nodes
#' arguments:
#' @param errorNode a extreme of the new edge
#' @param node another extreme of the new edge
#' @param edges edges defining the model
#' @param leftDeco deco for left node
#' @param rightDeco deco for right node
#' @return new list of edges after the addition
#'##############################################################
createEdge <- function(errorNode, node, edges, leftDeco, rightDeco){
  
  # gets the number of edges
  nrows <- nrow(edges)
  
  # adds one to nrow
  nrows <- nrows+1
  
  # adds the edge
  edges[nrows,] <- c(errorNode, node, leftDeco, rightDeco)
  
  # return edges
  edges
}

#'##############################################################
#' deletes a row from the data frame
#' arguments:
#' @param from first node
#' @param to second node
#' @param edges edges defining the model
#' @return new list of edges after removal
#'##############################################################
deleteEdge <- function(from, to, edges){
  
  # removes the edge
  index <- which((edges$from == from & edges$to == to) |
                   (edges$to == from & edges$from == to))
  edges <- edges[-index, ]
  
  # changes row names
  rownames(edges) <- 1:nrow(edges)
  
  # return edges
  return(edges)
}

#'#############################################################
#' deletes an edge selected at random between those belonging
#' to a given path
#' arguments:
#' @param path path to consider
#' @param edges set of edges to analyze
#' @return new list of edges after deletion
#'#############################################################
deleteEdgeAtRandom <- function(path, edges){
  # removes the cycle deleting an edge
  pathLength <- length(path)
  
  # selects a random number between 1 and pathLength
  index <- sample(1:pathLength,1)
  
  # deletes the edge
  from <- path[[index]]
  if (index < pathLength){
    to <- path[[index+1]]
  }
  else{
    to <- path[[1]]
  }
  
  # detele the edge
  edges <- deleteEdge(from, to, edges)
}

#'#############################################################
#' checks if there is an undirected edge between two nodes
#' arguments:
#' @param nodeA first node
#' @param nodeB second node
#' @param edges edges to consider
#' @return flag boolean value
#'#############################################################
checkUndirectedEdge <- function(nodeA, nodeB, edges){
  flag <- FALSE
  
  # gets undirected edges for both nodes
  undirectedEdges <- edges[((edges$from == nodeA & edges$to == nodeB &
                               edges$left == "none" & edges$right == "none") | 
                              (edges$from == nodeB & edges$to == nodeA &
                                 edges$left == "none" & edges$right == "none")),]
  
  if (!is.null(undirectedEdges) & nrow(undirectedEdges) != 0){
    flag <- TRUE
  }
  
  # return flag
  return(flag)
}

#'##############################################################
#' gets undirected edges for a given node
#' arguments:
#' @param node node to consider
#' @param edges edges where to look for
#' @return list of edges containing node
#'##############################################################
getUndirectedEdges <- function(node, edges){
  
  # gets undirected edges for node
  edgesForNode <- edges[((edges$from == node | edges$to == node) &
                           (edges$left == "none" & edges$right == "none")),]
}

#'#############################################################
#' get the nodes involved in a pair of edges, being node
#' the common node between them
#' arguments:
#' @param node common node between the edges
#' @param edge1 first edge
#' @param edge2 second edge
#' @return a list with nodes involved in both edges A - B - C
#'#############################################################
getNodesInEdges <- function(node, edge1, edge2){
  # sets nodeB
  nodeB <- node
  
  # get nodes A and C
  nodeA <- getOppositeNode(node, edge1)
  nodeC <- getOppositeNode(node, edge2)
  
  # return a list with these nodes
  return(list(A=nodeA, B=nodeB, C=nodeC))
}

#'##############################################################
#' get the opposite node to the one passed as first argument in the
#' edge passed as second argument
# arguments:
#' @param node reference node
#' @param edge edge to consider
#' @return node
#'##############################################################
getOppositeNode <- function(node, edge){
  other <- edge$from
  
  if (node == edge$from){
    other <- edge$to
  }
  
  # return node
  return(other)
}

#'##############################################################
#' get the spouses for a given node
#' arguments: 
#' @param node reference node
#' @param edges edges to analyze
#' @return list of nodes (spouses)
#'##############################################################
getSpousesForNode <- function(node, edges){
  # gets edges for node and with arraw-arrow as decoration
  edgesForNode <- edges[((edges$from == node | edges$to == node) & 
                           edges$left == "arrow" & edges$right == "arrow"), ]
  
  # consider every edge and insert the other node to spouses list
  spouses <- list()
  if (!is.null(edgesForNode) & nrow(edgesForNode) != 0){
    for(i in 1:nrow(edgesForNode)){
      edge <- edgesForNode[i,]
      
      # gets the other node
      spouses <- c(getOppositeNode(node, edge), spouses)
    }
  }
  
  # remove repetitions
  spouses <- unique(spouses)
  
  # return spouses
  return(spouses)
}

#'##############################################################
#' get spouses for every node
#' arguments:
#' @param edges edges to consider
#' @return list of spouses for each node
#'##############################################################
getSpouses <- function(edges){
  # gets all the nodes
  nodes <- unique(c(unique(edges$from), unique(edges$to)))
  
  # apply method for computing spouses for a given variable
  spouses <- lapply(nodes, getSpousesForNode, edges)
}

#'##############################################################
#' gets the neighnbours of node
#' arguments:
#' @param node node of reference
#' @param edges edges to analyze
#' @return list of neighbours
#'##############################################################
getNeighbours <- function(node, edges){
  # initializes neighbours
  neighbours <- list()
  
  # get the edges related to node
  edgesForNode <- edges[(edges$from == node | edges$to == node), ]
  
  # gets all the nodes in edges and removes node
  if (!is.null(edgesForNode) & nrow(edgesForNode) != 0){
    for(i in 1:nrow(edgesForNode)){
      neighbours <- c(getOppositeNode(node, edgesForNode[i,]), neighbours)
    }
  }
  
  # return neighbours
  return(neighbours)
}

#'##############################################################
#' gets the neighnbours of node but taking into account the direction.
#' if A -> B then B is neighbour of A, but A is not neighbour of B
#' arguments:
#' @param node reference node
#' @param edges edges to analyze
#' @return list of neighbours
#'##############################################################
getNeighboursWithDirections <- function(node, edges){
  
  # initializes neighbours
  neighbours <- list()
  
  # get the edges related to node
  edgesForNode <- edges[(edges$from == node | edges$to == node), ]
  
  # exclude edges with arrow in node side
  edgesForNode <- edgesForNode[!(((edgesForNode$to == node) & (edgesForNode$right == "arrow") &
                                    (edgesForNode$left == "none")) |
                                   ((edgesForNode$from == node) & (edgesForNode$left == "arrow") & 
                                      (edgesForNode$right == "none"))),]
  
  # gets all the nodes in edges and removes node
  if (!is.null(edgesForNode) & nrow(edgesForNode) != 0){
    for(i in 1:nrow(edgesForNode)){
      neighbours <- c(getOppositeNode(node, edgesForNode[i,]), neighbours)
    }
  }
  
  # return neighbours
  return(neighbours)
}

#'##############################################################
#' gets the neighnbours of node but taking into accounto only
#' undirected edges
#' arguments:
#' @param node reference node
#' @param edges edges to analyze
#' @return list of neighbours
#'##############################################################
getNeighboursWithUndirectedEdges <- function(node, edges){
  
  # initializes neighbours
  neighbours <- list()
  
  # get the edges related to node
  edgesForNode <- getUndirectedEdges(node, edges)
  
  # gets all the nodes in edges and removes node
  if (!is.null(edgesForNode) & nrow(edgesForNode) != 0){
    for(i in 1:nrow(edgesForNode)){
      neighbours <- c(getOppositeNode(node, edgesForNode[i,]), neighbours)
    }
  }
  
  # return neighbours
  return(neighbours)
}

#'##############################################################
#' method for getting the path between two nodes
#' arguments:
#' @param from start node
#' @param to destination node
#' @param visited flag of boolean values to control visited nodes
#' @param edges edges to analyze
#' @param neighboursFunction function to use for neighbours detection
#' @return list with three entries: boolean flag, list of visited
#'         nodes and list of non visited edges
#'##############################################################
getPath <- function(from, to=from, visited, edges, neighboursFunction){
  
  # inializes the result
  result <- list(flag=FALSE, path=visited, edges=edges)
  
  # base case: if to belongs to visited, return true
  if (any(visited == to)){
    result$flag <- TRUE
  }
  else{
    # inductive case: get neigbours of from
    neighbours <- neighboursFunction(from, edges)
    
    # gets the list of nodes to visit    
    toVisit <- setdiff(neighbours, visited)
    
    # consider every node to visit
    if(length(toVisit) != 0){
      for(i in 1:length(toVisit)){
        # sort nodes in lexicographical order
        toVisit <- sort(unlist(toVisit))
        
        # select the node
        nodeToVisit <- toVisit[i]
        
        # removes the edge between nodeToVisit and from
        edges <- deleteEdge(from, nodeToVisit, edges)
        
        # makes a recursive call
        result <- getPath(nodeToVisit, to, c(nodeToVisit, visited), edges, neighboursFunction)
        
        # if result is true, breaks the loop because the path
        # was found
        if (result$flag == TRUE){
          break
        }
      }
    }
  }
  
  # return result
  return(result)
}

#'####################### SECTION 2: MAMP functions #########################

#'##############################################################
#' check condition1 for MAMPCG
#' arguments: 
#' @param edges set of edges of the model
#' @return list with: boolean flag (true if the list of edges had
#'        to be changed) and list of resultant edges
#'##############################################################
checkCondition1 <- function(edges){
  changed <- FALSE
  continue <- TRUE
  
  while(continue){
    
    # get edges none - arrow
    candidateEdges <- edges[((edges$left == "none" & edges$right == "arrow") |
                               (edges$left == "arrow" & edges$right == "none")), ]
    
    # change continue value
    continue <- FALSE
    
    # checks the paths for every node
    if (nrow(candidateEdges) != 0){
      for(i in 1:nrow(candidateEdges)){
        # gets the edge
        edge <- candidateEdges[i,]
        
        # gets node from and to
        if (edge$left == "none"){
          nodeFrom <- edge$from
          nodeTo <- edge$to
        }
        else{
          nodeFrom <- edge$to
          nodeTo <- edge$from
        }
        
        result <- getPath(nodeTo, nodeFrom, list(nodeTo), edges, getNeighboursWithDirections)
        
        # if there is a path and the first edge is none - arrow
        # the condition 1 must be applied
        if (result$flag == TRUE){
          
          # removes an edge randomly selected
          edges <- deleteEdgeAtRandom(result$path, edges)
          
          # repeats the loop
          continue <- TRUE
          
          # sets changed to TRUE
          changed <- TRUE
          
          # breaks the for
          break
        }
      }
    }
  }
  return(list(changed=changed, edges=edges))
}

#'##############################################################
#' check condition2 for MAMPCG
#' arguments:
#' @param edges edges defining the model
#' @return list with: boolean flag (true if the list of edges had
#'        to be changed) and list of resultant edges
#'##############################################################
checkCondition2 <- function(edges){
  continue <- TRUE
  changed <- FALSE
  
  # detection loop
  while(continue){
    # get edges arrow - arrow
    candidateEdges <- edges[(edges$left == "arrow" & edges$right == "arrow"), ]
    
    # change continue value
    continue <- FALSE
    
    # checks the paths for every node
    if (nrow(candidateEdges) != 0){
      for(i in 1:nrow(candidateEdges)){
        # gets the edge
        edge <- candidateEdges[i,]
        
        # gets node from and to
        nodeFrom <- edge$from
        nodeTo <- edge$to
        
        # check if there is a path from nodeFrom to nodeTo with undirected edges
        result <- getPath(nodeTo, nodeFrom, list(nodeTo), edges, getNeighboursWithUndirectedEdges)
        
        # if there is a path and the first edge is none - arrow
        # the condition 1 must be applied
        if (result$flag == TRUE){
          
          # removes an edge randomly selected
          edges <- deleteEdgeAtRandom(result$path, edges)
          
          # repeats the loop
          continue <- TRUE
          
          # sets changed to TRUE
          changed <- TRUE
          
          # breaks the for
          break
        }
      }
    }
  }
  return(list(changed=changed, edges=edges))
}

#'##############################################################
#' check condition3 for MAMPCG
#' arguments:
#' @param edges edges defining the model
#' @return list with: boolean flag (true if the list of edges had
#'        to be changed) and list of resultant edges
#'##############################################################
checkCondition3 <- function(edges){
  changed <- FALSE
  continue <- TRUE
  
  # check loop
  while(continue){
    # sets continue to false
    continue <- FALSE
    
    # gets nodes
    nodes <- unique(c(edges$from, edges$to))
    
    # makes continue FALSE. Only with a change on the edges
    # this flag will be changed to TRUE
    continue <- FALSE 
    
    # considers every node
    for(i in 1:length(nodes)){
      node <- nodes[i]
      
      # gets undirected edges
      edgesForNode <- getUndirectedEdges(node, edges)
      
      # work only if there are at leas two edges
      if (!is.null(edgesForNode) & nrow(edgesForNode) >= 2){
        # considers every pair
        for(j in 1:(nrow(edgesForNode)-1)){
          for(k in (j+1):nrow(edgesForNode)){
            edge1 <- edgesForNode[j,]
            edge2 <- edgesForNode[k,]
            
            # gets the nodes involved in these edges: A - B - C
            nodesInEdges <- getNodesInEdges(node, edge1, edge2)
            
            # gets B node spouses
            bSpouses <- getSpousesForNode(node, edges)
            
            # if this set is not empty, then there must be an
            # endge between A and C
            if(!is.null(bSpouses) & length(bSpouses) != 0){
              # check the link between A and C
              flag <- checkUndirectedEdge(nodesInEdges$A, nodesInEdges$C, edges)
              
              # if flag is false, then adds an edge between A and C
              if (flag == FALSE){
                edges[(nrow(edges)+1),] <- c(from=nodesInEdges$A,to=nodesInEdges$C,
                                             left="none",right="none")
                
                # sets continue to TRUE
                continue <- TRUE
              }
            }
          }
        }
      }
    }
  }
  
  # return edges
  return(list(changed=changed, edges=edges))
}

#'##############################################################
#' check if the graph is a MAMPCG
#' arguments:
#' @param edges edges defining the model
#' @return list of edges required for a valid MAMPCG model
#'##############################################################
checkMAMPCG <- function(edges){
  
  # initializes flag to TRUE
  flag <- TRUE
  
  # while flag is TRUE
  while(flag){
    # check condition1
    res1 <- checkCondition1(edges)
    edges <- res1$edges
    
    # check condition2
    res2 <- checkCondition2(edges)
    edges <- res2$edges
    
    # check condition3    
    res3 <- checkCondition3(edges)
    edges <- res3$edges
    
    # compose the final result
    if (res1$changed == FALSE & res2$changed == FALSE & res3$changed == FALSE){
      flag=FALSE
    }
  }
  
  # return the set of edges
  return(edges)
}

#'##############################################################
#' generates a random graph with a certain probability for directed
#' undirected and bidirected graphs
#' arguments:
#' @param numberNodes number of nodes to consider
#' @param probDirected probability for directed edges
#' @param probUndirected probability for undirected edges
#' @param probBidirected probability for bidirected edges
#' @return list of resultant edges
#'##############################################################
generateRandomGraph <- function(numberNodes, probDirected, probUndirected, 
                                probBidirected){
  
  # generates a random graph
  rnet <- bnlearn::random.graph(LETTERS[1:numberNodes],method="ic-dag",
                                max.in.degree=2)
  
  # now gets the arcs
  rnetArcs <- bnlearn::arcs(rnet)
  
  # probability vector: probs for directed, undirected, bidirected
  probs <- c(probDirected, probUndirected, probBidirected)
  aprobs <- cumsum(probs)
  
  # generates a data frame with the required structure for edges
  edges <- data.frame(from=character(), to=character(), 
                      left=character(), right=character(), 
                      stringsAsFactors=FALSE)
  
  # considers every arc
  for(i in 1:nrow(rnetArcs)){
    # selects the edge
    arc <- rnetArcs[i,]
    
    # generates a random number
    rnumber <- runif(1)
    
    # gets the type according to rnumber
    type <- min(which(aprobs > rnumber))
    
    if (type == 1){
      # it is directed and nothing to do. Just insert the
      # edge
      edges[i,] <- c(arc["from"], arc["to"], "none","arrow")
    }
    else{
      if (type == 2){
        # it is undirected
        edges[i,] <- c(arc["from"], arc["to"], "none","none")
      }
      else{
        # bidirected
        edges[i,] <- c(arc["from"], arc["to"], "arrow","arrow")
      }
    }
  }
  
  # return edges
  return(edges)
}

#'##############################################################
#' method for generating a MAMPCG model
#' arguments:
#' @param numberNodes number of nodes
#' @param probs probs to use for the generation of directed, undirected
#'          and bidirected
#' @return list of edges defining the model
#'##############################################################
generateRandomMAMPCG <- function(numberNodes, probs){
  # generate the basic initial structure
  edges <- generateRandomGraph(numberNodes, probs[1], probs[2], probs[3])
  
  # checks the conditions
  edges <- checkMAMPCG(edges)
}

#'############### SECTION 3: functions for databases generation ################

#'##############################################################
#' method for transforming a set of edges in order to construct a
#' bayesian network
#' arguments:
#' @param edges edges defining the model
#' @return list with two entries: edges of the resultant BN and
#'         inmoralities produced by the conversion
#'##############################################################
transformToBayesianNetwork <- function(edges){
  # gets all the nodes
  nodes <- unique(c(edges$from,edges$to))
  
  # paste error prefix to every node
  rnodes <- sapply(nodes,function(node){
    rnode <- paste0("error",node)
  })
  
  # include an directed edge for errori to i
  for(i in 1:length(nodes)){
    # creates a new edge
    edges <- createEdge(rnodes[i],nodes[i],edges,"none","arrow")
  } 
  
  # selects undirected edges
  undirected <- edges[(edges$left == "none" & edges$right == "none"),]
  
  # initializes inmoralities to empty list
  inmoralities <- list()
  
  # remove the edge and add new edges between errorFrom and errorTo
  if (nrow(undirected) != 0){
    for(i in 1:nrow(undirected)){
      from <- undirected[i,]$from
      to <- undirected[i,]$to
      
      # removes the edge
      edges <- deleteEdge(from, to, edges)
      
      # add edges from error nodes to a new error node
      errorFrom <- paste0("error",from)
      errorTo <- paste0("error",to)
      error <- paste0("error",from)
      error <- paste0(error,to)
      edges <- createEdge(errorFrom, error, edges,"none","arrow")
      edges <- createEdge(errorTo, error, edges, "none", "arrow")
      
      # stores error into inmoralities list
      inmoralities <- unique(c(error, inmoralities))
    }
  }
  
  # selected bidirected edges
  bidirected <- edges[(edges$left == "arrow" & edges$right == "arrow"),]
  
  # for every bdirected node introduces links from common error
  # to error nodes
  if (nrow(bidirected) != 0){
    for(i in 1:nrow(bidirected)){
      from <- bidirected[i,]$from
      to <- bidirected[i,]$to
      
      # removes the edge
      edges <- deleteEdge(from, to, edges)
      
      # add edges from error to from and to
      error <- paste0("error",from)
      error <- paste0(error,to)
      edges <- createEdge(error, from, edges, "none", "arrow")
      edges <- createEdge(error, to, edges, "none", "arrow")
    }
  }
  
  # return edges and inmoralities
  return(list(edges=edges, inmoralities=inmoralities))
}

#'##############################################################
#' creates a Bnlearn net for helping the generation of distributions
#' arguments:
#' @param edges edges defining the model
#' @param check flag to show if the existance the cycles will be considered
#' @return resultant bnet
#'##############################################################
createBnlearnNet <- function(edges, check){
  # creates an empty graph with the variables
  nodes <- unique(c(edges$from, edges$to))
  
  # creates a bnet 
  bnet <- bnlearn::empty.graph(nodes)
  
  # now adds all the edges
  for(i in 1:nrow(edges)){
    edge <- edges[i,]
    #adds the arc
    bnet <- bnlearn::set.arc(bnet,from=edge$from,to=edge$to, check.cycles=check, 
                              debug=FALSE)
  }
  
  # return bnet
  return(bnet)
}

#'##############################################################
#' method for generating distributions for root nodes
#' arguments:
#' @param net net to be considered for databses generation
#' @return list of distributions for root nodes
#'##############################################################
generateDistributionForRootNodes <- function(net){
  # gets all the nodes without parents
  nodes <- bnlearn::nodes(net)
  
  # gets all the parents
  parentsOfNodes <- sapply(nodes,function(node){
    bnlearn::parents(net,node)
  })
  
  # initializes the list of distributions
  distributions <- list()
  
  # considers every node
  for (i in 1:length(nodes)){
    # gets node
    node <- nodes[i]
    
    # gets parents
    nodeParents <- parentsOfNodes[[i]]
    
    # check if there are no parents
    if (identical(nodeParents, character(0))){
      # gets a ramdom value between 1 and 2
      deviation <- runif(1)+1
      
      # sets the distribution
      distribution <- list(coef = c("(Intercept)" = 0), sd = deviation)
      
      # add the distribution to distributions
      distributions[[node]] <- distribution
    }
  }
  
  # return distributions
  distributions
}

#'##############################################################
#' method for generating distributions for non-root nodes
#' arguments:
#' @param net net to be considered for databses generation
#' @return list of distributions for non root nodes
#'##############################################################
generateDistributionForNonRootNodes <- function(net){
  # gets all the nodes without parents
  nodes <- bnlearn::nodes(net)
  
  # gets all the parents
  parentsOfNodes <- sapply(nodes,function(node){
    bnlearn::parents(net,node)
  })
  
  # initializes the list of distributions
  distributions <- list()
  
  # considers every node
  for (i in 1:length(nodes)){
    # gets node
    node <- nodes[i]
    
    # gets parents
    nodeParents <- parentsOfNodes[[i]]
    
    # node will have average = 0 and deviation = 0
    coefs <- c("(Intercept)"=0)
    
    # check if there are no parents
    if (!identical(nodeParents, character(0))){
      
      # considers every parent
      for(j in 1: length(nodeParents)){
        # generate the factor
        parent <- nodeParents[j]
        
        # checks if it is a error node
        if(length(grep("error",parent)) > 0){
          factor <- 1
        }
        else{
          factor <- runif(1)+1
        }
        
        # adds the factor to coefs
        coefs[parent] <- factor
      }
      
      # sets the distribution
      distribution <- list(coef=coefs, sd = 0)
      
      # add the distribution to distributions
      distributions[[node]] <- distribution
    }
  }
  
  # return distributions
  distributions
}

#'##############################################################
#' method for setting the parameters to the net
#' arguments: 
#' @param net net to consider
#' @param params parameters to set
#' @param resultant bnet.fit (with parameters)
#'##############################################################
setParameters <- function(net, params){
  # composes all the distributions
  net <- bnlearn::custom.fit(net, params)
}

#'##############################################################
#' method for remmoving from the complete data set the evidential
#' variables
#' arguments:
#' @param sample sample to filter removing evidence variables
#' @param evidenceVariables variables to remove from sample
#' @return dataframe with samples after removing columns for
#'         evidential variables
#'##############################################################
removeEvidentialVariables <- function(sample, evidenceVariables){
  # gets the complete list of nodes
  nodes <- names(sample)
  
  # removes evidence variables
  sample[ , -which(names(sample) %in% evidenceVariables)]
}

#'##############################################################
#' method for generating the complete dataset with the required sample
#' size
#' arguments: 
#' @param model model to sample from
#' @param sampleSize size of the sample to generate
#' @param threshold value to consider for evidence expressions
#' @param cl cluster to use (if possible to use several cores)
#' @return dataframe with samples
#'##############################################################
generateSample <- function(model, sampleSize, threshold, cl){
  # forms a expression where nodes in inmoralities are set to a value 
  # >= -threshold and <= threshold
  # forms the expressions for evidence: strB for >= expression
  # strL for <= expression and strC for the concatenation of both of them
  strB=paste("(",model$inmoralities, " >= ", -threshold, ")", sep="", collapse = " & ")
  strL=paste("(",model$inmoralities, " <= ", threshold, ")", sep="", collapse = " & ")
  strC=paste(strB, strL, sep=" & ")
  cat("Evidence expression: ", strC, "\n")

  # export strC to cluster nodes
  environment(strC) <- .GlobalEnv
  parallel::clusterExport(cl, "strC", envir=environment())
  
  # loop for getting the samples
  data <- NULL
  nSamples <- 0
  
  while(nSamples < sampleSize){
    # generate data. Perhaps the parameters must be changed for a faster generation
    # depending on the concrete model
    dataf <- bnlearn::cpdist(model$bnetSampling, nodes=bnlearn::nodes(model$bnet), 
                             evidence=eval(parse(text=strC)), method="ls", 
                             debug=FALSE, cluster=cl, batch=50000, 
                             n=sampleSize*30000)
 
    # updates the number of samples
    nSamples <- nSamples+nrow(dataf)
    cat("..",nSamples,"..")
    
    # join all the samples into data
    if(is.null(data)){
      data <- dataf
    }
    else{
      data <- rbind(data,dataf)
    }
  }

  # now remove all the variables containing error in their names
  data <- data[,-grep("error", colnames(data))]
  
  # remove extra samples
  data <- data[(1:sampleSize),]
  
  # return data
  return(data)
}

#'##############################################################
#' function for storing the data set to a file 
#' arguments:
#' @param sample: sample to store
#' @param id id of the database
#' @param numberSamples number of samples
#' @param path path where databases must be stored
#' @param filename filename to use
#' @return 
#'##############################################################
storeDatabase <- function(sample, id, numberSamples, path, filename){

  # compose the complete path to the file and creates the folder
  # if needed
  cpath <- paste0(path,filename)
  if (!file.exists(cpath)){
    dir.create(file.path(path,filename))
  }

  # gets the number of samples and concatenates with cpath
  if (!file.exists(paste(cpath,numberSamples,sep="/"))){
    dir.create(file.path(cpath,numberSamples))
  }
  cpath <- paste(cpath,numberSamples,sep="/")

  # now composes the name of the file
  filename <- paste(filename,numberSamples,sep="-")
  filename <- paste(filename,id,sep="-")
  filename <- paste0(filename,".db")
  filename <- paste(cpath,filename,sep="/")

  # writes the table
  write.table(sample,file=filename,col.names=TRUE, row.names=FALSE,sep=",")
}

#'##############################################################
#' method for storing the complete model: edges, bnet with
#' distributions and inmoralities
#' arguments:
#' @param model list containg all the information
#' @param folder where to store the model
#'##############################################################
storeModel <- function(model, folder){
  
  # compose the folder with model name
  pathName <- paste0(folder,model$name)
  pathName <- paste0(pathName,".mampcg")
  
  # uses the name of the model for storing info into an R file
  # (name, edges, bnet and inmoralities)
  saveRDS(model, pathName)
}

#'##############################################################
#' method for retrieving the complete model: edges, bnet with
#' distributions and inmoralities
#' arguments:
#' @param name of the model
#' @param folder where to look for the model
#' @return model
#'##############################################################
retrieveModel <- function(modelname, folder){
  
  # compose the folder with model name
  pathName <- paste0(folder,modelname)
  pathName <- paste0(pathName,".mampcg")
  
  # uses the name of the model for storing info into an R file
  # (name, edges, bnet and inmoralities)
  model <- readRDS(pathName)
}


#'##############################################################
#' prepare a MAMPCG model for sampling data from it
#' arguments:
#' @param edges
#' @param basename base name for the models
#' @param folder to analyze in order to assign an unique identifier
#' @return list with four entries: name, edges, bnet and inmoralities
#'##############################################################
prepareMAMPCGForSampling <- function(edges, basename, folder){
  
  # gets a unique name for this model
  name <- getNewFileName(basename, folder)
  
  # creates a BN from edges, without checking edges
  baseBnet <- createBnlearnNet(edges, check=FALSE)
  
  # transform into BN
  result <- transformToBayesianNetwork(edges)
  
  # create a bnlearn net for preparing parameters generaration
  bnet <- createBnlearnNet(result$edges, TRUE)
  
  # compute distributions for root nodes
  distributionsRootNodes <- generateDistributionForRootNodes(bnet)
  
  # compute distributions for non root nodes
  distributionsNonRootNodes <- generateDistributionForNonRootNodes(bnet)
  
  # compose the complete set of distributions
  distributions <- c(distributionsRootNodes, distributionsNonRootNodes)
  
  # now sets all the distributions to the net
  bnet <- setParameters(bnet, distributions)
  
  # finally return the name, edges, net and the set of inmoralities
  return(list(name=name, edges=edges, bnet=baseBnet, bnetSampling=bnet, inmoralities=result$inmoralities))
}

#'##############################################################
# method for generating data sets given a model
# arguments:
#' @param model model to use for generation
#' @param variants number of variants to generate
#' @param sampleSizes vector of sample sizes
#' @param threshold for evidence expressions
#' @param pathDb path where databases will be generated
#' @param cluster cluster to use several cores
#'##############################################################
generateDatabases <- function(model, variants, sampleSizes, threshold, pathDb, cluster){

  # considers every sample size
  for(ss in 1:length(sampleSizes)){
    
    # generates the variants
    for(v in 1:variants){
      # use this net for generating a sample
      data <- generateSample(model, sampleSizes[ss], threshold, cluster)
      
      # finally store the data
      storeDatabase(data, v, sampleSizes[ss], pathDb, model$name)
    }
  }
}

