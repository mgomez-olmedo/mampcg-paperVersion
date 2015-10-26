
library(bnlearn)
library(data.table)
library(R6)

#'##############################################################
#' class definition. Data members
#' moral: flag to define the way of considering v-structures
#' net: true Bayesian network
#' debug: flag for printing debug information
#' data: dataset used for learning
#' alpha: alpha for statistical tests of independence
#' pc: flag to adapt the algorithm in order to bahave as
#'     pc algorithm
#' originalEdges: true MAMP model
#' counters: counters for rules and steps application
#' edgesInfo: info about edges. This is needed to store
#' the required decoration for both sides of every edge
#' graph: graph with the structure of the model. Used for
#'        making easier the computation of adjacents
#' separators: list of separators between nodes
#'##############################################################
MAMPCGSearch <- R6Class(
  "MAMPCGSearch",
  public=list(
    #'##############################################################
    #' class contructor
    #' @param net
    #' @param data
    #' @param moral
    #' @param pc
    #' @param alpha
    #' @param edges
    #' @param debug
    #'############################################################## 
    initialize=function(net, data, moral=TRUE, pc=FALSE, alpha, edges=NULL, debug=FALSE){
      # sets debug flag
      private$debug <- debug
      
      # sets moral flag
      private$moral <- moral
      
      # sets net data member
      private$net <- net

      # stores the data
      private$data <- data
      
      # stores alpha value
      private$alpha <- alpha
      
      # flag showing the learn method will be normal
      # PC on dag
      private$pc <- pc
      
      # stores the info about the original structure defined
      # by edges and not by a net
      private$originalEdges <- edges

      # initializes the counters
      private$testCounter <- 0
      private$rule1Counter <- 0
      private$rule2Counter <- 0
      private$rule3Counter <- 0
      private$rule4Counter <- 0
      private$step10Counter <- 0
      private$step12Counter <- 0
      private$step13Counter <- 0
      private$step14Counter <- 0
      private$step15Counter <- 0
      
      # creates an empty data frame with columns named from, to, left,
      # and right
      private$edgesInfo <- data.frame(from=character(), to=character(), 
                                      left=character(), right=character(), 
                                      stringsAsFactors=FALSE)
      
      # creates an empty graph for these nodes
      private$graph <- bnlearn::empty.graph(names(data))
      
      # gets all the pairs of nodes
      combinations=combn(names(data), 2)
      
      # adds the set of undirected edges. The length of combinations is
      # divided by two just beacuse it contains pairs of elements. This is
      # achieved selecting all the columns
      for(i in 1:length(combinations[1,])){
        
        # adds the corresponding undirected arc
        private$edgesInfo[i,] <- c(combinations[1,i], combinations[2,i],"none","none")
        
        # add the edge
        private$graph <- bnlearn::set.edge(private$graph,from=combinations[1,i],to=combinations[2,i],
                                  check.cycles=FALSE, debug=FALSE)
      }
      
      # creates an empty list for separators
      private$separators <- list()
    },
    
    #'##############################################################
    #' compares structures between this object and the true
    #' net (stored in private$net)
    #' moral argument to show if v-structures must consider triplexes
    #' A -> B <- C with A and C as adjacents or not (TRUE this v-structure
    #' is considered; with false this is not a valid v-structure)
    #'##############################################################
    compareStructure = function(){
      if (private$debug){
        cat("\n -------------------------- compareStructure begin ---------------------\n")
      }
      
      # compares the structures
      # gets the skeleton of the true net
      sk <- bnlearn::empty.graph(nodes(private$net))
      bnlearn::amat(sk,ignore.cycles=TRUE) <- bnlearn::amat(private$net)
      trueDependenceGraph <- bnlearn::skeleton(sk)

      # gets th skeleton of the current graph
      sk <- bnlearn::empty.graph(nodes(private$graph))
      bnlearn::amat(sk,ignore.cycles=TRUE) <- bnlearn::amat(private$graph)
      currentDependenceGraph <- bnlearn::skeleton(sk)

      # se realiza la comparacion
      result <- compare(trueDependenceGraph,currentDependenceGraph,arcs=FALSE)
      
      # now computes the precision and recall measurements
      recall <- result$tp/(result$tp+result$fn)
      precision <- result$tp/(result$tp+result$fp)
      cat("Recall: ",recall," precision: ",precision,"\n")

      # now get the vstructures of the true net
      if (is.null(private$originalEdges)){
        trueVsdf <- as.data.frame(vstructs(as.bn(modelstring(private$net)), moral=private$moral), 
                                  stringsAsFactors=FALSE)
        
        # adds a new column with the type
        trueVsdf[,"type"] <- c(rep(1,nrow(trueVsdf)))
        if (private$debug){
          cat("V-structures in true net (computed with net information): \n")
          print(trueVsdf)
          cat("---------------------------------------------\n")
        }
      }
      else{
        trueVsdf <- private$getModelVStructuresDataFrame(private$originalEdges)
        if (private$debug){
          cat("V-structures in original net (computed with edges information): \n")
          print(trueVsdf)
          cat("---------------------------------------------\n")
        }
      }
            
      # gets th vstructures of the learnt model
      currentVsdf <- private$getModelVStructuresDataFrame(private$edgesInfo)
      if (private$debug){
        cat("V-structures in learnt net: \n")
        print(currentVsdf)
        cat("---------------------------------------------\n")
      }
      
      # now compares both sets of vstructures
      resultVs <- private$compareVStructures(trueVsdf,currentVsdf)
      recallVs <- resultVs$recall
      precisionVs <- resultVs$precision
      cat("RecallVs: ",recallVs," precisionVs: ",precisionVs,"\n\n")
      
      if (private$debug){
        cat("\n -------------------------- compareStructure end ---------------------\n")
      }
      
      # check for NAs and change with 0
      if (is.na(precision)){
         precision=0
      }
      if (is.na(precisionVs)){
         precisionVs=0
      }
      
      # finally returns all this measurements
      return(list(recall=recall, precision=precision, recallVs=recallVs, precisionVs=precisionVs))
    },
    
    #'##############################################################
    #' sets edgesInfo information. This is used only for testing
    #' @param edges
    #'##############################################################
    setEdgesInfo = function(edges){
      private$edgesInfo <- edges
    },
    
    #'##############################################################
    #' sets graph data member. This is used only for testing
    #' @param graph
    #'##############################################################
    setGraph = function(graph){
      private$graph <- graph
    },
    
    #'##############################################################
    #' tests rule1 application. Used for testing
    #'##############################################################
    testRule1 = function(){
      if (private$debug){
        cat("\n -------------------------- testRule1 begin ---------------------\n")
      }
      
      # applies initialLoop
      private$initialLoop()
      
      if (private$debug){
        cat("Initial loop end: apply rule 1 phase ................... \n")
      }
      
      # show info about bInfo with separators info
      self$show(TRUE)
      
      # applies rule 1
      flag <- TRUE
      while(flag){
        flag <- private$applyRule1()
      }
      
      if (private$debug){
        cat("\n -------------------------- testRule1 end ---------------------\n")
      }
    },
    
    #'##############################################################
    #' tests rule2 application. Used for testing
    #'##############################################################
    testRule2 = function(){
      if (private$debug){
        cat("\n -------------------------- testRule2 begin ---------------------\n")
      }
      
      # show info about bInfo with separators info
      self$show(TRUE)
      
      # applies rule 1
      flag <- TRUE
      while(flag){
        flag <- private$applyRule2()
      }
      
      if (private$debug){
        cat("\n -------------------------- testRule2 end ---------------------\n")
      }
    },
    
    #'##############################################################
    #' tests rule3 application. Used for testing
    #'##############################################################
    testRule3 = function(){
      private$applyRule3()
    },
    
    #'##############################################################
    #' tests rule 4 application. Used for testing
    #'##############################################################
    testRule4 = function(){
      private$applyRule4()
    },
    
    #'##############################################################
    #' tests the application of the step 10
    #'##############################################################
    testStep10 = function(){
      # applies initialLoop
      private$initialLoop()
      
      result <- TRUE
      while(result == TRUE){
        result <- private$changeEdgesInCycle()
      }
    },
    
    #'##############################################################
    #' tests the learn algorithm without making independence tests. Only
    #' for testing purposes
    #'##############################################################
    testLearn = function(){
      if (private$debug){
        cat("Begin of rules 1 - 2 - 3 - 4 application\n")
      }
            
      # loop of rules application
      rules <- c(TRUE,TRUE,TRUE,TRUE)
      private$applyRules(rules)
      
      
      if (private$debug){
        cat("End of rules 1 - 2 - 3 - 4 application\n")
      }
      
      if (private$debug){
        cat("Step 10 application\n")
      }
      
      # applies step 10
      result <- TRUE
      while(result){
        result <- private$changeEdgesInCycle()
      }
      
      if (private$debug){
        cat("End of step 10 application\n")
      }
      
      if (private$debug){
        cat("Begin of rules 2 - 3 - 4 application\n")
      }
      
      # applies rules 2, 3 and 4
      # loop of rules application
      rules <- c(FALSE,TRUE,TRUE,TRUE)
      private$applyRules(rules)
      
      if (private$debug){
        cat("End of rules 2 - 3 - 4 application\n")
      }
      
      if (private$debug){
        cat("Step 12 application\n")
      }
      
      # now changes edges block-wildcard by none-arrow (step 12)
      counter <- private$changeEdgesDecoration("block","none","none","arrow")
      private$step12Counter <- private$step12Counter+counter
      
      if (private$debug){
        cat("End of step 12 application\n")
      }
      
      if (private$debug){
        cat("Step 13 application\n")
      }
      
      # now changes edges block-block by arrow-arrow (step 13)
      counter <- private$changeEdgesDecoration("block","block","arrow","arrow")
      private$step13Counter <- private$step13Counter+counter
      
      if (private$debug){
        cat("End of tep 13 application\n")
      }
      
      if (private$debug){
        cat("Step 14 application\n")
      }
      
      # apply step 14 
      result <- TRUE
      while(result){
        result <- private$applyStep14()
      }
      
      if (private$debug){
        cat("End of step 14 application\n")
      }
      
      if (private$debug){
        cat("Step 15 application\n")
      }
      
      # apply step 15 while true
      result <- TRUE
      while(result){
        result <- private$applyStep15()
      }
      
      if (private$debug){
        cat("End of step 15 application\n")
      }
    },
    
    #'##############################################################
    #' tests the decoration changes. Only for testing
    #'##############################################################
    testChangeEdgesDecoration = function(){
      private$changeEdgesDecoration("block","block","none","none")
    },
    
    #'##############################################################
    #' add separators info just for testing. Only for testing
    #' @param pair pair of nodes
    #' @param separatorSet separators for nodes in pair
    #'##############################################################
    addSeparator = function(pair, separatorSet){
      private$storeSeparatorInfo(pair,separatorSet)
    },
    
    #'##############################################################
    #' main method for learning
    #'##############################################################
    learn=function(){
      if (private$debug){
        cat("\n -------------------------- learn begin ---------------------\n")
        cat("PC mode: ",private$pc,"\n")
      }
      
      # applies initialLoop
      private$initialLoop()
            
      if (private$debug){
        cat("Begin of rules 1 - 2 - 3 - 4 application\n")
      }
      
      # loop of rules application
      rules <- c(TRUE,TRUE,TRUE,TRUE)
      private$applyRules(rules)
      
      
      if (private$debug){
        cat("End of rules 1 - 2 - 3 - 4 application\n")
      }
      
      # step 10 only works in mampc mode
      if (!private$pc){
        if (private$debug){
          cat("Step 10 application\n")
        }
      
        # applies step 10
        result <- TRUE
        while(result){
          result <- private$changeEdgesInCycle()
        }
      
        if (private$debug){
          cat("End of step 10 application\n")
        }
      }
      
      # step 11 only works in mampc mode
      if (!private$pc){
        if (private$debug){
          cat("Begin of rules 2 - 3 - 4 application\n")
        }
      
        # applies rules 2, 3 and 4
        # loop of rules application
        rules <- c(FALSE,TRUE,TRUE,TRUE)
        private$applyRules(rules)
      
        if (private$debug){
          cat("End of rules 2 - 3 - 4 application\n")
        }
      }
      
      if (private$debug){
        cat("Step 12 application\n")
      }
      
      # now changes edges block-wildcard by none-arrow (step 12)
      counter <- private$changeEdgesDecoration("block","none","none","arrow")
      private$step12Counter <- private$step12Counter+counter
      
      if (private$debug){
        cat("End of step 12 application\n")
      }
      
      # step 13 only works in mampc mode
      if (!private$pc){
        if (private$debug){
          cat("Step 13 application\n")
        }
      
        # now changes edges block-block by arrow-arrow (step 13)
        counter <- private$changeEdgesDecoration("block","block","arrow","arrow")
        private$step13Counter <- private$step13Counter+counter
      
        if (private$debug){
          cat("End of tep 13 application\n")
        }
      }
      
      # step 14 only works in mampc mode
      if (!private$pc){
        if (private$debug){
          cat("Step 14 application\n")
        }
      
        # apply step 14 
        result <- TRUE
        while(result){
          result <- private$applyStep14()
        }
      
        if (private$debug){
          cat("End of step 14 application\n")
        }
      }
      
      # step 14 only works in mampc mode
      if (!private$pc){
        if (private$debug){
          cat("Step 15 application\n")
        }
      
        # apply step 15 while true
        result <- TRUE
        while(result){
          result <- private$applyStep15()
        }
      
        if (private$debug){
          cat("End of step 15 application\n")
        }
      }
      
      # finally takes teh changes to graph in pc mode
      if(private$pc){
        private$updateGraph()
      }
      
      if (private$debug){
        cat("\n -------------------------- learn end ---------------------\n")
      }
    },
        
    #'##############################################################
    #' method for showing the data contained by the object
    #'##############################################################
    show = function(showSeparators=FALSE){
      if (private$debug){
        cat("\n -------------------------- BuildInfo object ---------------------\n")
        cat("PC mode: ",private$pc,"\n")
        cat("Edges info: \n")
        cat("Number of edges: ",nrow(private$edgesInfo),"\n")
        print(private$edgesInfo)
        cat("....................................................................\n")
        cat("Associated bn: \n")
        print(private$graph)
        cat("Nodes in graph: ",nodes(private$graph),"\n")
        cat("Show separators info: ",showSeparators,"\n")
        if (showSeparators == TRUE){
          cat("....................................................................\n")
          cat("Separators set: \n")
          print(private$separators)
          cat("....................................................................\n")
        }
        cat("Tests counter: ",private$testCounter,"\n")
        cat("Rule 1 counter: ",private$rule1Counter,"\n")
        cat("Rule 2 counter: ",private$rule2Counter,"\n")
        cat("Rule 3 counter: ",private$rule3Counter,"\n")
        cat("Rule 4 counter: ",private$rule4Counter,"\n")
        cat("Step 10 counter: ",private$step10Counter,"\n")
        cat("Step 12 counter: ",private$step12Counter,"\n")
        cat("Step 13 counter: ",private$step13Counter,"\n")
        cat("Step 14 counter: ",private$step14Counter,"\n")
        cat("Step 15 counter: ",private$step15Counter,"\n")
      }
    }
  ),
  
  #'##############################################################
  #' PRIVATE SECTION OF THE CLASS
  #'##############################################################  
  # all the data are kept private as well as some auxiliar functions
  private = list(
    net="bn.fit",
    data="data.frame",
    moral="logical",
    alpha="numeric",
    pc="logical",
    edgesInfo="data.frame",
    graph="bn",
    separators="list",
    originalEdges="data.frame",
    debug="logical",
    testCounter="numeric",
    rule1Counter="numeric",
    rule2Counter="numeric",
    rule3Counter="numeric",
    rule4Counter="numeric",
    step10Counter="numeric",
    step12Counter="numeric",
    step13Counter="numeric",
    step14Counter="numeric",
    step15Counter="numeric",
    
    #'##############################################################
    #' method performing the first stage of the algorithm
    #'##############################################################
    initialLoop=function(){
      if (private$debug){
        cat("\n -------------------------- initialLoop begin ---------------------\n")
        self$show()
      }
      
      # initializes l to 0
      l <- 0
      
      # gets the number of nodes
      nnodes <- length(nodes(private$graph))
      
      # main loop
      while(l <= (nnodes-2)){
        if (private$debug){
          cat("Loop for l: ",l,"\n")
        }
        
        # forms the pairs of nodes
        pairs <- combn(bnlearn::nodes(private$graph),2)
        # forms the pairs of nodes
        pairs <- combn(bnlearn::nodes(private$graph),2)
        pairs1 <- lapply(seq_len(ncol(pairs)), function(i) {c(pairs[1,i],pairs[2,i])})
        pairs2 <- lapply(seq_len(ncol(pairs)), function(i) {c(pairs[2,i],pairs[1,i])})
        pairs <- append(pairs1,pairs2)
        
        # gets the pairs with the proper conditions
        filteredPairs <- lapply(pairs, private$treatPair, l=l)
        
        # adds 1 to l
        l <- l+1
      }
      if (private$debug){
        cat("\n -------------------------- initialLoop end ---------------------\n")
      }
    },
    
    #'##############################################################
    #' private method for making the operations with a pair of nodes to detect
    #' independence/dependence
    #'##############################################################
    treatPair=function(pair, l){
      if (private$debug){
        cat("\n -------------------------- treatPair begin ---------------------\n")
        cat("Pair to be tested for l =",l, " \n")
        print(pair)
      }
      
      # check the conditions in order to apply the test
      result <- private$checkConditionsTestAnalysis(pair, l)
      
      # if the flag of result is TRUE, it is needed to make
      # independence tests
      if (result$flag == TRUE){
        # make the tests of this pair
        private$makeTests(pair,l,result$set)
      }
      
      if (private$debug){
        cat("\n -------------------------- treatPair end ---------------------\n")
      }
    },
        
    #'##############################################################
    #' function for making the independence tests on a pair of nodes
    #'##############################################################
    makeTests=function(pair,l,set){
      if (private$debug){
        cat("\n -------------------------- makeTests begin ---------------------\n")
      }
      
      # makes all the sets of size equals to l
      sets <- combn(set,l)
      
      # if l is 0, sets will be be a list with an empty list
      if (length(sets) == 0){
        sets <- list(list())
      }
      else{
        sets <- lapply(seq_len(ncol(sets)), function(i) {
          separatorSet <- list()
          for(j in 1:l){
            separatorSet[j] <- sets[j,i]
          }
          separatorSet
        })
      }
      
      # now considers every set and make the corresponding independence
      # test. If the rest result is positive, change the graph and the
      # info about edges. This must be done in a loop style in order to
      # avoid useless tests
      if (private$debug){
        cat("Set of possible separators: \n")
        for(i in 1:length(sets)){
          cat(" ",unlist(sets[[i]]), " ")
        }
        cat("\n")
        cat("Number of candidates: ",length(sets),"\n")
      }
      
      for(i in 1:length(sets)){
        # gets the separator
        separator <- as.character(sets[[i]])
        
        #if (length(separator) == 0){
        #  separator <- NULL
        #}
        
        if (private$debug){
          cat("Longth of separator: ", length(separator),"\n")
          cat("Separator: |",separator,"|\n")
          cat("x: ")
          print(pair[[1]])
          cat("y: ")
          print(pair[[2]])
        }
        
        # makes the test
        if (!is.null(private$originalEdges)){
          # continuous case
          testResult <- bnlearn::ci.test(x=pair[[1]], y=pair[[2]],z=separator,private$data,test="mi-g")
        }
        else{
          # discrete case
          testResult <- bnlearn::ci.test(x=pair[[1]], y=pair[[2]],z=separator,private$data,test="mi")
        }
        
        # increments the counter of test
        private$testCounter <- private$testCounter+1
        
        if (private$debug){
          cat("Test for ",pair[[1]]," - ",pair[[2]], " separator: ")
          print(separator)
          cat("test result: \n")
          print(testResult)
        }
        
        # if p-value >= alpha, considers independence
        if (testResult$p.value >= private$alpha){
          if (private$debug){
            cat("removing edge...............\n")
          }
          
          # removes the edge
          private$removeEdge(pair)
          
          # stores information about separator
          if(is.null(separator)){
            separator <- list()
          }
          private$storeSeparatorInfo(pair,separator)
          
          # in order to debug the code, shows the result
          if(private$debug){
            cat("--------------------- Current state --------------------\n")
            self$show()
            cat("--------------------------------------------------------\n")
            cat("\n")
          }
          
          # breaks the loop
          break
        }
      }
      
      if (private$debug){
        cat("\n -------------------------- makeTests end ---------------------\n")
      }
    },
    
    #'##############################################################
    #' method to remove an edge in the graph
    #' arguments:
    #' @param pair pair of nodes defining the edge. It is not possible
    #'         the existance of several edges between a pair of nodes
    #'##############################################################
    removeEdge=function(pair){
      if (private$debug){
        cat("\n -------------------------- removeEdge begin ---------------------\n")
        cat("Removing edge ",pair[[1]], " - ",pair)
      }
      
      # removes the edge from the graph
      private$graph <- bnlearn::drop.edge(private$graph,pair[[1]],pair[[2]])
      
      # removes the info about this edge
      private$edgesInfo <- private$edgesInfo[!(private$edgesInfo$from == pair[[1]] &
                                                 private$edgesInfo$to == pair[[2]]),]
      
      if (private$debug){
        cat("\n -------------------------- removeEdge end ---------------------\n")
      }
    },
    
    #'##############################################################
    #' method to store information about separator
    #' arguments: 
    #' @param pair: pair of nodes under analysis
    #' @param separtor: separator set for nodes in pair
    #'##############################################################
    storeSeparatorInfo=function(pair, separator){
      if (private$debug){
        cat("\n -------------------------- storeSeparator begin ---------------------\n")
        cat("Separator for ",pair[[1]], " - ",pair[[2]],"\n")
        cat("Separator set: \n")
        print(separator)
      }
      
      # forms both keys
      keyab=paste(pair[[1]],pair[[2]],sep="")
      keyba=paste(pair[[2]],pair[[1]],sep="")
      
      # stores the info
      private$separators[[keyab]] <- separator
      private$separators[[keyba]] <- separator
      if (private$debug){
        cat("\n -------------------------- storeSeparator end ---------------------\n")
      }
    },
    
    #'##############################################################
    #' method for applying the rules once and again while possible
    #' arguments:
    #' @param rules logical vector with boolean flags for rule1, rule2,
    #'         rule3 and rule4
    #'##############################################################
    applyRules = function(rules){
      if (private$debug){
        cat("\n -------------------------- applyRules begin ---------------------\n")
      }
      keep <- TRUE
      # Sets flags to false to avoid problems if the rules will not
      # be applied
      flag1 <- FALSE
      flag2 <- FALSE
      flag3 <- FALSE
      flag4 <- FALSE
      
      # apply every rule as much as possible
      while(keep){
        if (private$debug){
          cat("....................... rule application loop ......................\n")
          self$show(TRUE)
        }
        
        # applies rule 1
        if (rules[1] == TRUE){
          flag1 <- private$applyRule1()
        }
        
        # applies rule 2
        if (rules[2] == TRUE){
          flag2 <- private$applyRule2()
        }
        
        # applies rule 3
        if (rules[3] == TRUE){
          flag3 <- private$applyRule3()
        }
        
        # applies rule 3
        if (rules[4] == TRUE){
          flag4 <- private$applyRule4()
        }
        
        # updates keep value
        keep <- flag1 | flag2 | flag3 | flag4
        
        if (private$debug){
          cat("flag1: ",flag1, " flag2: ", flag2, " flag3: ", flag3, " flag4: ",flag4,"\n")
          cat("Value of keep in rules application: ",keep,"\n")
        }
      }
      
      if (private$debug){
        cat("\n -------------------------- applyRules end ---------------------\n")
        cat("flag1: ",flag1, " flag2: ", flag2, " flag3: ", flag3, " flag4: ",flag4,"\n")
        cat("------------------------------------------------------------------\n")
      }
    },
    
    #'##############################################################
    #' method for rule 1 application
    #'##############################################################
    applyRule1 = function(){
      if (private$debug){
        cat("\n -------------------------- applyRule1 begin ---------------------\n")
      }

      # flag to show if the rule was applied
      applied <- FALSE
      
      # get pairs with a common node
      pairs <- private$matchEdgesWithCommonNode(private$edgesInfo)
      
      # filter pairs by separators: in this case B must not be included in
      # the separator set for A and B
      pairs <- private$filterPairsBySeparators(pairs, FALSE)
            
      # # generate a data frame if there are pairs; in any other case
      # NULL will be returned
      candidatePairsdf <- private$generateAlignmentDataFrame(pairs)
      
      # checks if there a pair and selects pairs of edges with the required pattern
      if (private$containsData(candidatePairsdf)){
        # selects the pattern block - wildcard - wildcard - wildcard
        candidatePairsdf <- private$selectPairPattern(candidatePairsdf, "wildcard","wildcard", 
                                                      "wildcard", "noblock")
                
        # considers every candidate and breaks the loop if one is found
        if (private$containsData(candidatePairsdf)){
          for(i in 1:nrow(candidatePairsdf)){
            # selects the corresponding pair of edges
            pairInfo <- candidatePairsdf[i, ]
            
            # check the existence of an directed edge between the extremes
            acEdge <- private$checkEdgeInArgumentEdges(pairInfo$anode,pairInfo$cnode,private$edgesInfo)
                   
            # change decoration and breaks the loop
            if (!acEdge){
              # change edge info: first edge becomes block-wildcard and second
              # becomes wildcard-block.
              changed1 <- private$changeEdgeDecoration(pairInfo$anode, pairInfo$bnode, 
                                           "block", "wildcard")
              changed2 <- private$changeEdgeDecoration(pairInfo$bnode, pairInfo$cnode,
                                           "wildcard", "block")
              
              if (changed1 & changed2){
                # adds 1 to rule1Counter
                private$rule1Counter <- private$rule1Counter+1
              
                # sets applied be true
                applied <- TRUE
              
                # breaks the loop
                break
              }
              else{
                # restore the edge changed
                if (changed1){
                  private$changeEdgeDecoration(pairInfo$anode,pairInfo$bnode,
                                                 pairInfo$adeco,pairInfo$bdeco1)
                }
                if (changed2){
                  private$changeEdgeDecoration(pairInfo$bnode,pairInfo$cnode,
                                                 pairInfo$bdeco2,pairInfo$cdeco)
                }
              }
            }
          }
        }
      }
      
      if (private$debug){
        cat("\n -------------------------- applyRule1 end ---------------------\n")
        cat("Applied: ",applied,"\n")
        cat("\n ---------------------------------------------------------------\n")
      }
      
      # return applied
      return(applied)
    },
    
    #'##############################################################
    #' method for rule 2 application
    #'##############################################################
    applyRule2 = function(){
      if (private$debug){
        cat("\n -------------------------- applyRule2 begin ---------------------\n")
      }
      
      # flag to show if the rule was applied
      applied <- FALSE
      
      # get pairs with a common node
      pairs <- private$matchEdgesWithCommonNode(private$edgesInfo)
      
      # filter pairs by separators: in this case B must be included in
      # the separator set for A and B
      pairs <- private$filterPairsBySeparators(pairs, TRUE)
      
      # # generate a data frame if there are pairs; in any other case
      # NULL will be returned
      candidatePairsdf <- private$generateAlignmentDataFrame(pairs)
      
      # checks if there a pair
      if (private$containsData(candidatePairsdf)){
        # selectes the pattern block - wildcard - wildcard - wildcard
        candidatePairsdf <- private$selectPairPattern(candidatePairsdf, "block",
                                              "wildcard", "noblock", "wildcard")
        
        # considers every candidate and breaks the loop if one if found
        if (private$containsData(candidatePairsdf)){
          for(i in 1:nrow(candidatePairsdf)){
            pairInfo <- candidatePairsdf[i, ]
            
            # check the existence of an directed edge between the extremes
            acEdge <- private$checkEdgeInArgumentEdges(pairInfo$anode,pairInfo$cnode,private$edgesInfo)
            
            # if flag is true change edges decoration on the extremes
            # not in common
            if (!acEdge){
              # change edge info: first edge remains and the second becomes to
              # block-wildcard
              changed <- private$changeEdgeDecoration(pairInfo$bnode, pairInfo$cnode,
                                   "block", "wildcard")
              
              if (changed){
                # adds 1 to rule2Counter
                private$rule2Counter <- private$rule2Counter+1
              
                # sets 1 to applied
                applied <- TRUE
              
                # break the loop
                break
              }
            }
          }
        }
      }
      
      if (private$debug){
        cat("\n -------------------------- applyRule2 end ---------------------\n")
        cat("Applied: ",applied,"\n")
        cat("\n ---------------------------------------------------------------\n")
      }
      
      # return applied
      return(applied)
    },
    
    #'##############################################################
    #' method for rule 3 application. The change in edges and the increment
    #' in the counter are made in auxiliar methods
    #'##############################################################
    applyRule3 = function(){
      if (private$debug){
        cat("\n -------------------------- applyRule3 begin ---------------------\n")
      }
      # considers from-to direction
      flag <- private$applyRule3WithDirection(FALSE)
      
      # if the rule was not applied, considers the opposite direction
      if (flag == FALSE){
        flag <- private$applyRule3WithDirection(TRUE)
      }
      
      if (private$debug){
        cat("\n -------------------------- applyRule3 end ---------------------\n")
      }
      
      # return flag
      return(flag)
    },
    
    #'##############################################################
    #' method for rule 4 application
    #'##############################################################
    applyRule4=function(){
      if (private$debug){
        cat("\n -------------------------- applyRule4 begin ---------------------\n")
      }
      
      # begins with the flag set to FALSE
      applied <- FALSE
      
      # gets the pairs with a common node
      pairs <- private$matchEdgesWithCommonNode(private$edgesInfo)
      
      # generates a data frame if there are pairs; in any other case
      # NULL will be returned
      candidatePairsdf <- private$generateAlignmentDataFrame(pairs)
      
      # checks if there are some pair
      if(private$containsData(candidatePairsdf)){
        # selects the target pattern: wildcard-wildcard-block-wildcard
        candidatePairsdf <- private$selectPairPattern(candidatePairsdf, "wildcard","wildcard",
                                                      "block","wildcard")
        
        if (private$debug){
          cat("Pair of edges with wildcard - wildcard - block - wildcard pattern: ",nrow(candidatePairsdf),"...............\n")
          print(candidatePairsdf)
        }
        
        # now get sets of two pairs in order to get the rombus
        # pattern for thr rule 4 application
        if (private$containsData(candidatePairsdf) & nrow(candidatePairsdf) > 1){
          # gets the quartets
          quartets <- private$getQuartets(candidatePairsdf)
          
          # checks if quartets is null
          if (private$containsData(quartets)){
            if (private$debug){
              cat("Quartets for rule application: ",nrow(quartets),"...............\n")
              print(quartets)
            }
            
            # if there are quartets, check if the rule can be applied,
            # being a in the separator set of CD
            if (length(quartets) != 0){
              quartet <- quartets[[1]][[1]]
              result <- private$checkConditionForRule4(quartet)
              
              # if flag is true then the link between quartet$pair1$anode and
              # quartet$pair2$cnode can ne changed
              if (result$flag){
                if (private$debug){
                  print(result)
                }
                
                # changes the edge
                changed <- private$changeEdgeDecoration(quartet$pair1$anode, quartet$pair1$cnode,"block","wildcard")
              
                if (changed){
                  # sets applied to true
                  applied <- TRUE
                
                  # adds 1 to rule4Counter
                  private$rule4Counter <- private$rule4Counter+1
                }
              }
            }
          }
        }
      }
      
      if (private$debug){
        cat("\n -------------------------- applyRule4 end ---------------------\n")
      }
      
      # return applied
      return(applied)
    },
    
    #'##############################################################
    #' method for the application of step 10
    #'##############################################################
    changeEdgesInCycle = function(){
      if (private$debug){
        cat("\n -------------------------- changeEdgesInCycle begin ---------------------\n")
      }
      
      # initialize the process
      lnodes <- bnlearn::nodes(private$graph)
      nnodes <- length(lnodes)
      flag <- FALSE
      
      # initially the list of edges is the whole set
      edges  <- private$edgesInfo  
      
      # it is needed to consider every node
      for(i in 1:nnodes){
        
        # selects the node
        node <- lnodes[[i]]
        
        if (private$debug){
          cat("Looking for loops from node ",node,"\n")
        }
        
        # checks if the node belongs to a cycle
        result <- private$isReachableWithLoopAndDeco(node,node,"noblock","noblock",
                                                     list(), edges)
        
        if (private$debug){
          cat("Result on investigate cycles from ",node,"\n")
          cat("result: ",result$flag,"\n")
          cat("Path: \n")
          print(result$path)
          cat("Length of path: ",result$length,"\n")
        }
        
        # is is reachable in a path without blocks and with length
        # bigger or equal to 4, it is needed
        # to change all the edges to block-block
        if (result$flag == TRUE){   
          # first at all checks the induced subgraph
          subgraph <- private$getInducedSubgraph(result$path)
          if (private$debug){
            cat("Path with loop detected: \n")
            print(subgraph)
          }
          
          if (result$length > 3 & nrow(subgraph) == result$length){
            flag <- TRUE
            
            # change edges in path
            private$changeEdgesInPath(result$path,"block","block") 
            
            # adds 1 to the counter
            private$step10Counter <- private$step10Counter+1
            
            # break the loop
            break  
          }
          else{              
            # remove the edges contained in path
            edges <- private$removeEdgesInArgumentInPath(result$path, edges)
            
            if (private$debug){
              cat("Reduced the set of edges for further exploration: \n")
              print(edges)
            }
          }            
        }
        else{
          # if there is a loop, remove the edges in it
          if (result$loop == TRUE){
            if (private$debug){
              cat("Loop detection for previous nodes\n")
              print(result$path)
            }
            edges <- private$removeEdgesInArgumentInCycle(result$path, edges)
          }
        }
      }
      
      if (private$debug){
        cat("\n -------------------------- changeEdgeInCycle end ---------------------\n")
        cat("Flag value to return: ",flag,"\n")
      }
      
      # return the flag
      return(flag)
    },
    
    #'##############################################################
    #' method for changing the decorations of the edges following a
    #' certain pattern: the two first decorations are changed by
    #' the new ones (needed for steps 12 and 13)
    #' arguments:
    #' @param deco1 actual decoration of one side
    #' @param deco2 actual decoration of other side
    #' @param newDeco1 new decoration for side1
    #' @param newDeco2 new decoration for side2
    #'##############################################################
    changeEdgesDecoration = function(deco1, deco2, newDeco1, newDeco2){
      if (private$debug){
        cat("\n -------------------------- changeEdgesDecoration begin ---------------------\n")
        cat("deco1: ",deco1," deco2: ",deco2,
            " newDeco1: ",newDeco1, " newdDeco2: ",newDeco2,"\n")
      }
      
      # sets counter to 0
      counter <- 0
      
      # get the edges with deco1 in left and deco2 in right
      directedEdges <- private$getEdgesWithPatternFromTo(deco1,deco2)

      # changes the edges decoration
      if (private$containsData(directedEdges)){
        counter <- counter+nrow(directedEdges)
        # changes all of them
        apply(directedEdges, 1, function(edge){
          edge <- unlist(edge)          
          private$changeEdgeDecoration(edge["from"],edge["to"],newDeco1,newDeco2)
        })
        
        if (private$debug){
          cat("inverted edges in from - to direction\n")
        }
      }
            
      # the same with the set of inverted edges
      invertedEdges <- private$getEdgesWithPatternFromTo(deco2,deco1)
      
      # changes the edges decoration
      if (private$containsData(invertedEdges)){
        counter <- counter+nrow(invertedEdges)
        # change all of them
        apply(invertedEdges, 1, function(edge){
          edge <- unlist(edge)        
          private$changeEdgeDecoration(edge["from"],edge["to"],newDeco2,newDeco1)
        })
        
        if (private$debug){
          cat("inverted edges in to - from direction\n")
        }
        
      }
      
      if (private$debug){
        cat("\n -------------------------- changeEdgesDecoration end ---------------------\n")
      }
      
      # return flag
      return(counter)
    },
    
    #'##############################################################
    #' method for step14 application
    #'##############################################################
    applyStep14 = function(){
      if (private$debug){
        cat("\n -------------------------- applyStep14 begin ---------------------\n")
      }
      
      # flag to show if the rule was applied
      applied <- FALSE
      
      # get pairs of edges with common node
      pairs <- private$matchEdgesWithCommonNode(private$edgesInfo)
      
      if (private$debug){
        cat("Pairs of edges with common node: \n")
        print(pairs)
      }
            
      # considers the pairs
      if (private$containsData(pairs)){
                
        # selects pairs with pattern arrow - arrow - none - none
        candidatePairsdf <- private$selectPairPattern(pairs, "arrow", "arrow",
                                                      "arrow", "arrow")
        
        if (private$debug){
          cat("Candidate pairs with pattern arrow - arrow - arrow - arrow \n")
          print(candidatePairsdf)
        }
        
        # considers every pair
        if (private$containsData(candidatePairsdf)){
          for(i in 1:nrow(candidatePairsdf)){
            pairInfo <- candidatePairsdf[i,]
            if (private$debug){
              cat("Candidate pair for step 14 application: \n")
              print(pairInfo)
            }
          
            # checks the existence of a directed edge between the extremes
            # in the induced subgraph 
            acEdge <- private$checkEdgeInArgumentEdges(pairInfo$anode,pairInfo$cnode,private$edgesInfo)
          
            # change decoration and breaks the loop
            if (!acEdge){
              # change edge info: both edges become none - none
              private$changeEdgeDecoration(pairInfo$anode, pairInfo$bnode, 
                                         "none", "none")
              private$changeEdgeDecoration(pairInfo$bnode, pairInfo$cnode, 
                                         "none", "none")
              
              print("Se suma 1 a step14Counter......")
            
              # adds 1 to rule1Counter
              private$step14Counter <- private$step14Counter+1
            }
          }
        }
      }
      
      if (private$debug){
        cat("\n -------------------------- applyStep14 end ---------------------\n")
        cat("Applied: ",applied,"\n")
        cat("\n ---------------------------------------------------------------\n")
      }
      
      # return applied
      return(applied)
    },
    
    #'##############################################################
    #' method for step15 application
    #'##############################################################
    applyStep15 = function(){
      if (private$debug){
        cat("\n -------------------------- applyStep15 begin ---------------------\n")
      }
      
      # flag to show if the rule was applied
      applied <- FALSE
      
      # get pairs with a common node
      pairs <- private$matchEdgesWithCommonNode(private$edgesInfo)
      
      # # generate a data frame if there are pairs; in any other case
      # NULL will be returned
      candidatePairsdf <- private$generateAlignmentDataFrame(pairs)
      
      # checks if there a pair
      if (private$containsData(candidatePairsdf)){
        # selectes the pattern block - wildcard - wildcard - wildcard
        candidatePairsdf <- private$selectPairPattern(candidatePairsdf, 
                                                      "arrow", "arrow", "none", "none")
        
        # considers every candidate and breaks the loop if one if found
        if (private$containsData(candidatePairsdf)){
            for(i in 1:nrow(candidatePairsdf)){
              pairInfo <- candidatePairsdf[i, ]
            
              # this method also performs the changes in the edges
              # in order to prepare the change in edgesInfo
              directedEdge <- private$checkEdgeForNodesWithPattern(pairInfo$anode,pairInfo$cnode,"none","none")
            
              # if flag is true change edges decoration on the extremes
              # not in common
              if (directedEdge){
                # change the second edge
                private$changeEdgeDecoration(pairInfo$anode,pairInfo$bnode,"none","none")
              
                # set applies to true
                applied <- TRUE
              
                # adds 1 to the counter
                private$step15Counter <- private$step15Counter+1
              
                # break the loop
                break
             }
           }
        }
      }
      
      if (private$debug){
        cat("\n -------------------------- applyStep15 end ---------------------\n")
        cat("Applied: ",applied,"\n")
        cat("\n ---------------------------------------------------------------\n")
      }
      
      # return applied
      return(applied)
    },
    
    #'##############################################################
    #' method for changing an edge between two nodes. The new
    #' edge decorations are stored in the argument
    #' arguments:
    #' @param edge edge eith information about involved nodes and
    #'         new decorations
    #'##############################################################
    changeEdge = function(edge){
      if (private$debug){
        cat("\n -------------------------- changeEdge begin ---------------------\n")
        cat("Changing  edge:\n")
        print(edge)
      }
      
      # looks for the index of the corresponding edge in edgesInfo
      index <- which(private$edgesInfo$from == edge$from & private$edgesInfo$to == edge$to)
      cat("Indice: ",index)
      
      # if it is found, change it
      if(length(index) == 1){
        private$edgesInfo[index,] <- edge
      }
      
      if (private$debug){
        cat("\n -------------------------- changeEdge end ---------------------\n")
      }
    },
    
    #'##############################################################
    #' method for changing the edges int the path
    #' arguments:
    #' @param path path where changes will be applied
    #' @param deco1 decoration of a side
    #' @param deco2 decoration of the other side
    #'##############################################################
    changeEdgesInPath = function(path, deco1, deco2){
      if (private$debug){
        cat("\n -------------------------- changeEdgesInPath begin ---------------------\n")
        cat("Path: \n")
        print(path)
        cat("deco1: ",deco1," deco2: ",deco2)
      }
      
      # considers every edge
      for(i in 1:length(path)){
        index1 <- i
        index2 <- i+1
        if (index2 > length(path)){
          index2=1
        }
        
        # gets the next two nodes
        anode <- path[[index1]]
        bnode <- path[[index2]]
        
        # change the decorations to block-block
        private$changeEdgeDecoration(anode,bnode,deco1,deco2)
      }
      if (private$debug){
        cat("\n -------------------------- changeEdgesInPath end ---------------------\n")
      }
    },
        
    #'##############################################################
    #' look for edges containing node (as from or to) and having  
    #' nodeDecorator on its side and otherDecorator in the opposite
    #' one
    #' arguments:
    #' @param node node to look for in edges
    #' @param nodeDecorator the method will work on edges having node
    #                  and this decoration on a side
    #' @param otherDecorator decoration for the other side
    #' @param edges list of edges to analyze
    #'##############################################################
    getEdgesInArgumentForNode = function(node, nodeDecorator, otherDecorator, edges){
      if (private$debug){
        cat("\n -------------------------- getEdgesInArgumentForNode begin ---------------------\n")
        cat("Node: ",node, "  node decorator: ",nodeDecorator, " other decorator: ",otherDecorator," \n")
        cat("List of edges: \n")
        print(edges)
      }
      
      # get edges outgoing from node
      outgoingEdges <- edges[edges$from == node,]
      
      # considers left decoration if is not wildcard
      if (nodeDecorator == "noblock"){
        outgoingEdges <- outgoingEdges[outgoingEdges$left != "block",]
      } 
      else{
        if (nodeDecorator != "wildcard"){
          outgoingEdges <- outgoingEdges[outgoingEdges$left == nodeDecorator,]
        }
      }
      
      # considers right decoration
      if (otherDecorator == "noblock"){
        outgoingEdges <- outgoingEdges[outgoingEdges$right != "block",]
      }
      else{
        if (otherDecorator != "wildcard"){
          outgoingEdges <- outgoingEdges[outgoingEdges$right == otherDecorator,]
        }
      }
      
      if (private$debug){
        cat("List of outgoing edges detected: \n")
        print(outgoingEdges)
      }
      
      # get edges reaching node
      incomingEdges <- edges[edges$to == node,]
      
      # considers node decoration
      if (nodeDecorator == "noblock"){
        incomingEdges <- incomingEdges[incomingEdges$right != "block",]
      }
      else{
        if (nodeDecorator != "wildcard"){
          incomingEdges <- incomingEdges[incomingEdges$right == nodeDecorator,]
        }
      }
      
      # considers otherNode decoration
      if (otherDecorator == "noblock"){
        incomingEdges <- incomingEdges[incomingEdges$left != "block",]
      }
      else{
        if (otherDecorator != "wildcard"){
          incomingEdges <- incomingEdges[incomingEdges$left == otherDecorator,]
        }
      }
      
      if (private$debug){
        cat("List of incoming edges detected: \n")
        print(incomingEdges)
      }
      
      # compound the final list
      result <- merge(outgoingEdges, incomingEdges, all=TRUE)
      
      if (private$debug){
        cat("\n -------------------------- getEdgesInArgumentForNode end ---------------------\n")
        cat("Complete list of edges: \n")
        print(result)
        cat("-----------------------------------------------------------------------\n")
      }
      
      # return result
      return(result)
    },
    
    #'##############################################################
    #' matches edges with a common node: align the edges and form a
    #' new data structure to be used for rules applications
    #' arguments:
    #' @param edges list of edges to consider
    #'##############################################################
    matchEdgesWithCommonNode = function(edges){
      # make cartesian product
      prod <- merge(edges,edges,all=TRUE,by=NULL)
      
      # remove matches made with the same edge
      prodf <- prod[!(prod$from.x == prod$from.y & prod$to.x == prod$to.y),]
      
      # select pairs with a common node
      # e1 - e2 (to-from)
      # e1 - e2(R) (to-to)
      # e1(R) - e2 (from - from)
      # e1(R) - e2(R) (from - to)
      prodfSelected <- prodf[(prodf$to.x == prodf$from.y |
                                prodf$to.x == prodf$to.y |
                                prodf$from.x == prodf$from.y |
                                prodf$from.x == prodf$to.y),]
      
      # return the result
      return(prodfSelected)
    },
    
    #'##############################################################
    #' generate a data frame from the list of lists with alignment pairs. This
    #' is an auxiliar method for rules application
    #' arguments: 
    #' @param pairs pairs of nodes to consider (list of lists)
    #'##############################################################
    generateAlignmentDataFrame = function(pairs){
      if (private$debug){
        cat("\n -------------------------- generateAlignmentDataFrame begin ---------------------\n")
        cat("Number of pairs: ",nrow(pairs),"\n")
      }
      
      result <- NULL
      # checks if there is any pair of edges with a common node

      if (private$containsData(pairs)){
         # call align edges on pairs
         index <- seq(1:nrow(pairs))
        
         # generate the list with the result
         candidatePairs <- lapply(index, private$alignEdges, df=pairs)
         result <-data.table::rbindlist(candidatePairs)
      }
      
      if (private$debug){
        cat("\n -------------------------- generateAlignmentDataFrame end ---------------------\n")
        cat("Final candidate pairs: \n")
        if (!is.null(result)){
          print(result)
        }
        else{
          cat("no pairs with required conditions.....\n")
        }
      }
      
      # return result
      return(result)
    },
    
    #'##############################################################
    #' method for aligning the edges with a common node between them
    #' in order to get groups of nodes candidates for rules application
    #' This is an auxiliar method for rules application
    #' arguments:
    #' @param index index of the pair of edges to consider
    #' @param df data structure with the aligned pairs of edges
    #'##############################################################
    alignEdges = function(index, df){  
      row <- df[index,]
      
      # get the nodes in both edges
      nodes1 <- c(row$from.x, row$to.x)
      nodes2 <- c(row$from.y, row$to.y)
      
      # gets B (coomon node)
      bNode <- intersect(nodes1, nodes2)
      
      # gets the difference
      diff <- setdiff(union(nodes1,nodes2),bNode)
      aNode <- diff[1]
      cNode <- diff[2]
      
      # assigns the decorations: it is needed to check the
      # pattern of the alignment
      # e1 - e2 (to-from)
      # e1 - e2(R) (to-to)
      # e1(R) - e2 (from - from)
      # e1(R) - e2(R) (from - to)
      e1Reverted <- FALSE
      e2Reverted <- FALSE
      if (row$to.x == bNode){
        # first edge os ok
        adeco <- row$left.x
        bdeco1 <- row$right.x
      }
      else{
        # first edge is reverted
        adeco <- row$right.x
        bdeco1 <- row$left.x
        e1Reverted <- TRUE
      }
      
      # check for the second edge
      if (row$from.y == bNode){
        # compounds the sequence of decorations
        bdeco2 <- row$left.y
        cdeco <- row$right.y
      }
      else{
        # the second edge is reverted
        bdeco2 <- row$right.y
        cdeco <- row$left.y
        e2Reverted <- TRUE
      }
      
      # compound the info
      return(list(anode=aNode, bnode=bNode, cnode=cNode, 
                  adeco=adeco, bdeco1=bdeco1, bdeco2=bdeco2, cdeco=cdeco,
                  e1reverted=e1Reverted, e2reverted=e2Reverted))
    },
    
    #'##############################################################
    #' from a data frame with aligned edges, selects those matching a certain
    #' pattern adeco- bdeco1- bdeco2 - cdeco
    #' arguments:
    #' @param pairs of aligned edges
    #' @param deco1
    #' @param deco2
    #' @param deco3
    #' @param deco4
    #'##############################################################
    selectPairPattern = function(pairs, deco1, deco2, deco3, deco4){
      if (private$debug){
        cat("\n -------------------------- selectPairPattern begin ---------------------\n")
        cat("Decorations to look for: ",deco1, " - ",deco2, " -", deco3, " - ", deco4, "\n")
      }
      
      # initializes selected with pairs
      selected <- pairs 
      
      # considers deco1
      if (deco1 == "noblock"){
        selected <- selected[selected$adeco != "block",]
      }
      else{
        if (deco1 != "wildcard"){
          selected <- selected[selected$adeco == deco1,]
        }
      }
      
      # considers deco2
      if (deco2 == "noblock"){
        selected <- selected[selected$bdeco1 != "block",]
      }
      else{
        if (deco2 != "wildcard"){
          selected <- selected[selected$bdeco1 == deco2,]
        }
      }
      
      # considers deco3
      if (deco3 == "noblock"){
        selected <- selected[selected$bdeco2 != "block",]
      }
      else{
        if (deco3 != "wildcard"){
          selected <- selected[selected$bdeco2 == deco3,]
        }
      }
      
      # considers deco4
      if (deco4 == "noblock"){
        selected <- selected[selected$cdeco != "block",]
      }
      else{
        if (deco4 != "wildcard"){
          selected <- selected[selected$cdeco == deco4,]
        }
      }
      
      if (private$debug){
        cat("Selected pairs:  \n")
        print(selected)
        cat("\n -------------------------- selectPairPattern end ---------------------\n")
      }
      
      # return selected
      return(selected)
    },
    
    #'##############################################################
    #' check the conditions for a pair of nodes in order to apply independence 
    #' tests on it
    #' arguments:
    #' @param pair pair of nodes
    #' @param l size of separator set to consider
    #'##############################################################
    checkConditionsTestAnalysis = function(pair, l){
      if (private$debug){
        cat("\n -------------------------- checkConditionTestAnalysis begin ---------------------\n")
      }
      
      # initializes result
      result <- list(flag=FALSE)
      
      # gets a and b nodes
      anode <- pair[[1]]
      bnode <- pair[[2]]
      
      # gets a adjacents
      aAdjacents <- private$getAdjacents(anode)
      
      # gets b adjacents
      bAdjacents <- private$getAdjacents(bnode)
      
      # check if anode belongs to bAdjacents
      if (anode %in% bAdjacents){
        
        # the set of adjacents to consider will depend on the pc mode
        # in mampc mode gets the adjacents of the adjacents set. In PC
        # mode only considers aAdjacents
        if (!private$pc){
          # gets ad(ad(anodes))
          aAdjacentsAdjacents <- private$getAdjacentsSet(aAdjacents)
        }
        else{
          aAdjacentsAdjacents <- aAdjacents
        }
        
        # removes a and b
        aAdjacentsAdjacents <- aAdjacentsAdjacents[aAdjacentsAdjacents != anode &
                                                     aAdjacentsAdjacents != bnode]
        
        # now checks the lenght of aAdjacents
        if (length(aAdjacentsAdjacents) >= l){
          result <- list(pair=list(anode,bnode), flag=TRUE, set=aAdjacentsAdjacents)
        }
      }
      
      if (private$debug){
        cat("\n -------------------------- checkConditionTestAnalysis end ---------------------\n")
        cat("Pair: ",pair[[1]]," - ", pair[[2]], " flag: ",result$flag, 
            " candidate set: ",result$set,"\n")
        cat("----------------------------------------------------------------------------------\n")
      }
      
      # return result
      return(result)
    },
        
    #'##############################################################
    #' method for computing the adjacents for a node
    #' arguments: 
    #' @param node node to consider
    #'##############################################################
    getAdjacents = function(node){
      # return the neighbour of the node
      return(bnlearn::nbr(private$graph,node))
    },
    
    #'##############################################################
    #' method for computing the adjacents of a set of nodes
    #' finally append adjacents to the computed set
    #' @param set os nodes to get adjacents 
    #'##############################################################
    getAdjacentsSet = function(set){  
      # get the adjacents of all the nodes in setsesion2.zip
      setAdjacents <- unlist(lapply(set,function(x){private$getAdjacents(x)}))
      setAdjacents <- unique(append(set,setAdjacents))
    },
    
    #'##############################################################
    #' check if there is a directed edge between the extreme nodes
    #' arguments:
    #' @param pairInfo info about the pair of edges to consider
    #'##############################################################
    checkDirectedEdgeBetweenExtremes = function(pairInfo){
      # initializes variables
      directedEdge <- FALSE
      
      if (private$debug){
        cat("\n -------------------------- checkDirectedEdgeBetweenExtremes begin ---------------------\n")
        cat("A: ",pairInfo$anode, " B: ",pairInfo$bnode, " C: ", pairInfo$cnode,"\n")
        cat("Pattern: ",pairInfo$adeco, " - ", pairInfo$bdeco1, " - ", 
            pairInfo$bdeco2, " - ", pairInfo$cdeco,"\n")
      }
      
      # removes the edges with other nodes
      involvedNodes <- list(pairInfo$anode, pairInfo$bnode, pairInfo$cnode)
      inducedEdges <- private$getInducedSubgraph(involvedNodes)
      
      # checks if there is an edge between anode and cnode in inducedEdges
      directedEdge <- private$checkEdgeInArgumentEdges(pairInfo$anode, pairInfo$cnode, inducedEdges)
      if (private$debug){
        cat("Direct edges between ",pairInfo$anode, " and ",pairInfo$cnode," : ",directedEdge,"\n")
      }
      
      if (private$debug){
        cat("\n -------------------------- checkDirectedEdgeBetweenExtremes end ---------------------\n")
        cat("Directed edge: ", directedEdge,"\n")
        cat("------------------------------------------------------------------------------\n")
      }
      
      # return flag
      return(directedEdge)
    },
    
    #'##############################################################
    #' check if a pair of edges has the conditions for applying rule 2
    #' arguments:
    #' @param pairInfo information about the pair of edges to analyze
    #'##############################################################
    checkConditionForRule2 = function(pairInfo){
      flag <- FALSE
      
      if (private$debug){
        cat("\n -------------------------- checkConditionForRule2 begin ---------------------\n")
        cat("A: ",pairInfo$anode, " B: ",pairInfo$bnode, " C: ", pairInfo$cnode,"\n")
        cat("Pattern: ",pairInfo$adeco, " - ", pairInfo$bdeco1, " - ", 
            pairInfo$bdeco2, " - ", pairInfo$cdeco,"\n")
      }
      
      # removes the edges with other nodes
      involvedNodes <- list(pairInfo$anode, pairInfo$bnode, pairInfo$cnode)
      inducedEdges <- private$getInducedSubgraph(involvedNodes)
      
      # checks if there is an edge between anode and cnode in inducedEdges
      directEdge <- private$checkEdgeInArgumentEdges(pairInfo$anode, pairInfo$cnode, inducedEdges)
      if (private$debug){
        cat("Direct edges between ",pairInfo$anode, " and ",pairInfo$cnode," : ",directEdge,"\n")
      }
      
      # flag is just the opposite to directEdge
      flag <- !directEdge
            
      if (private$debug){
        cat("\n -------------------------- checkConditionForRule2 end ---------------------\n")
      }
      
      # return flag
      return(flag)
    },
    
    #'##############################################################
    #' method for applying rule3 with a given direction given by inverted argument
    #' arguments:
    #' @param inverted boolean flag with the folllowing meaning
    #'           inverted - false: direction from - to
    #'           inverted - true: direction to - from
    #'##############################################################
    applyRule3WithDirection=function(inverted){
      if (private$debug){
        cat("\n -------------------------- applyRule3WithDirection begin ---------------------\n")
      }
      
      # flag to show if the rule was applied
      applied <- FALSE
      direction <- NULL
      
      # if direction is from-to then looks for noblock-wildcard edges
      # in to-from direction looks for wildcard-noblock
      if (inverted == FALSE){
        # gets noblock-wildcard edges
        candidateEdges <- private$edgesInfo[private$edgesInfo$left != "block",]
        if (private$debug){
          # print the list of candidate edges
          direction <- "Direction from - to\n"
          cat(direction)
          print(candidateEdges)
        }
      }
      else{
        # gets wildcard-noblock
        candidateEdges <- private$edgesInfo[private$edgesInfo$right != "block",]
        
        if (private$debug){
          # print the list of candidate edges
          direction <- "Direction to - from\n"
          cat(direction)
          print(candidateEdges)
        }
      }
            
      # before considering a candidate edge between A-B it is needed
      # to check that: a) there is at least an edge with block in A
      # side and b) there is at least an edge with block in B opposite
      # side. This tries to avoid useless checks
      if (private$containsData(candidateEdges)){
        flags <- apply(candidateEdges,1,private$checkInitialConditionForRule3, inverted=inverted)
        candidateEdges <- candidateEdges[flags,]
        #candidateEdges <- private$checkInitialConditionForRule3(candidateEdges, inverted)
      }
      
      if (private$debug){
        cat("Rule 3 candidates (",direction,"): ",nrow(candidateEdges),"...............\n")
        print(candidateEdges)
      }
      
      # proceeds to consider every candidate edge (between A and B)
      if (private$containsData(candidateEdges)){
        # apply rule 3 on these candidate edges
        applied <- private$applyRule3OnCandidateEdges(candidateEdges, inverted)
      }
      
      if (private$debug){
        cat("\n -------------------------- applyRule3WithDirection end ---------------------\n")
        cat("Applied: ",applied, "\n")
        cat("------------------------------------------------------------------\n")
      }
      
      # return applied
      return(applied)
    },
    
    #'##############################################################
    #' apply rule 3 on candidate edges
    #' arguments:
    #' @param candidateEdges edges to consider for applying rule 3
    #' @param direction: boolean flag showing the direction to consider
    #'              in the edges
    #'##############################################################
    applyRule3OnCandidateEdges=function(candidateEdges, direction){
      if (private$debug){
        cat("\n -------------------------- applyRule3OnCandidateEdges begin ---------------------\n")
        cat("Direction: ", direction)
      }
      
      # flag to show if the rule was applied
      applied <- FALSE
      
      # considers every pair
      for(i in 1:nrow(candidateEdges)){
        candidateEdge <- candidateEdges[i,]
        
        if(private$debug){
          cat("Candidate edge for rule 3 application: \n")
          print(candidateEdge)
          cat("......................................................\n")
        }
        
        # selects origin and dest according to the direction
        if (direction == FALSE){
          origin <- candidateEdge$from
          dest <- candidateEdge$to
        }
        else{
          origin <- candidateEdge$to
          dest <- candidateEdge$from
        }
        
        # checks if candidateEdge$to is reachable from candidateEdge$from
        # with block - wildcard. Block must be on from side and wildcard
        # on destination side (to). The edges are passed as argument int
        # order to hava a recursive behaviour
        result <- private$isReachableWithLoopAndDeco(origin, dest, "block", "wildcard", 
                                           list(origin), private$edgesInfo)
        
        if(private$debug){
          cat("Node ",candidateEdge$to, " reachable from ", 
              candidateEdge$from,": ",result$flag,"\n") 
        }
        
        # if the flag is true, then change the arc and break the loop
        if (result$flag){
          changed <- private$changeEdgeDecoration(origin, dest, "block", "wildcard")
          
          if (changed){
            # updates applied
            applied <- TRUE
          
            # adds 1 to rule3Counter
            private$rule3Counter <- private$rule3Counter+1
          
            if (private$debug){
              cat("............ Rule 3 applied ............\n")
            }
          
            break
          }
        }
      }
      
      if (private$debug){
        cat("\n -------------------------- applyRule3OnCandidateEdges end ---------------------\n")
        cat("Applied: ",applied, "\n")
        cat("------------------------------------------------------------------\n")
      }
      
      # return applied
      return(applied)
    },
    
    #'##############################################################
    #' checks if the edge A-B has additional edges in order
    #' to be a candidate for rule 3 application
    #' arguments:
    #' @param edge edge to consider
    #' @param inverted boolean flag 
    #'##############################################################
    checkInitialConditionForRule3 = function(edge, inverted){
      if (private$debug){
        cat("\n -------------------------- checkInitialConditionForRule3 begin ---------------------\n")
        cat("Edge under analysis: ")
        print(edge)
        cat("Inverted: ",inverted,"\n")
      }
      
      flag <- FALSE
      
      # checks if the edge is a candidate for rule3 application:
      if (inverted == FALSE){
        # there must be an edge with block in from side
        # there must be an edge with block in the opposite side of to side
        edges1 <- private$edgesInfo[(private$edgesInfo$from == edge["from"] &
                                       private$edgesInfo$left == "block") |
                              (private$edgesInfo$to == edge["from"] &
                                 private$edgesInfo$right == "block"),]
        
        edges2 <- private$edgesInfo[(private$edgesInfo$to == edge["to"] &
                                       private$edgesInfo$left == "block") |
                              (private$edgesInfo$from == edge["to"] &
                                 private$edgesInfo$right == "block"),]
        
        if (private$debug){
          cat("Number of edges with block on A ",edge["from"]," side: ",nrow(edges1),"\n")
          cat("Number of edges with block on B ",edge["to"]," opposite side: ",nrow(edges2),"\n")
        }
        
        # gives value to flag
        if (private$containsData(edges1) & private$containsData(edges2)){
          flag <- TRUE
        }
      }
      else{
        # checks the opposite direction
        edges1 <- private$edgesInfo[(private$edgesInfo$from == edge["to"] &
                                       private$edgesInfo$left == "block") |
                             (private$edgesInfo$to == edge["to"] &
                                private$edgesInfo$right == "block"),]
        edges2 <- private$edgesInfo[(private$edgesInfo$to == edge["from"] &
                                       private$edgesInfo$left == "block") |
                             (private$edgesInfo$from == edge["from"] &
                                private$edgesInfo$right == "block"),]
        if (private$debug){
          cat("Number of edges with block on A ",edge["from"]," side: ",nrow(edges1),"\n")
          cat("Number of edges with block on B ",edge["to"]," opposite side: ",nrow(edges2),"\n")
        }
        
        if (!is.null(edges1) & !is.null(edges2) & nrow(edges1) != 0 & nrow(edges2) != 0){
          flag <- TRUE
        }
      }
      
      if (private$debug){
        cat("\n -------------------------- checkInitialConditionForRule3 end ---------------------\n")
        cat("Flag returned: ",flag,"\n")
      }
      
      # return flag
      return(flag)
    },
    
    #'##############################################################
    #' check if a pair of edges has the conditions for applying rule 4
    #' arguments:
    #' @param quartet quartet of edges candidates to rule 4 application
    #'##############################################################
    checkConditionForRule4 = function(quartet){
      flag <- FALSE
      edge1 <- quartet$pair1
      edge2 <- quartet$pair2
      
      if (private$debug){
        cat("\n -------------------------- checkConditionForRule4 begin ---------------------\n")
        cat("A: ",edge1$anode, " C: ",edge1$bnode, " B: ", edge1$cnode,"\n")
        cat("Pattern: ",edge1$adeco, " - ", edge1$bdeco1, " - ", 
            edge1$bdeco2, " - ", edge1$cdeco,"\n")
        cat("A: ",edge2$anode, " D: ",edge2$bnode, " B: ", edge2$cnode,"\n")
        cat("Pattern: ",edge2$adeco, " - ", edge2$bdeco1, " - ", 
            edge2$bdeco2, " - ", edge2$cdeco,"\n")
      }
      
      # there must be a link noblock - wildcard between edge1$anode and
      # edge1@cnode (entre A y B)
      directLink <- private$checkEdgeWithDecorations(edge1$anode, edge1$cnode, "noblock","wildcard")
      
      # if there is such a link, keeps on checking the condition
      if (directLink){
        
        # first at all get the induces subgraph
        involvedNodes <- list(edge1$anode,edge1$bnode,edge1$cnode,edge2$bnode)
        inducedEdges <- private$getInducedSubgraph(involvedNodes)
        
        # checks the presence of tha direct link between C and D (ege1$bnode and
        # edge2$bnode)
        cdEdge <- private$checkEdgeInArgumentEdges(edge1$bnode,edge2$bnode,private$edgesInfo)
        
        # keeps on testing if there are no such link
        if (!cdEdge){
          
          # now check if B belongs to AC separator
          # b should not belongs to ac separators
          cdSeparator <- private$separators[[paste(edge1$bnode,edge2$bnode,sep="")]]
        
          if (private$debug){
            cat("Separator set between ",edge1$bnode," and ", edge2$bnode,"\n")
            print(cdSeparator)
            cat("Edge to change: \n")
            print(directLink)
          }
        
          # If everything is ok return the flag
          if (edge1$anode %in% cdSeparator){
            flag <- TRUE
          }
        }
      }
      if (private$debug){
        cat("\n -------------------------- checkConditionForRule4 end ---------------------\n")
      }
      
      # return flag
      return(list(flag=flag,edge=directLink))
    },
    
    #'##############################################################
    #' check if a pair of edges has the conditions for applying step14
    #' arguments:
    #' @param pairInfo info about the pair of edges to consider
    #'##############################################################
    checkConditionForStep14 = function(pairInfo){
      flag <- FALSE
      
      if (private$debug){
        cat("\n -------------------------- checkConditionForStep14 begin ---------------------\n")
        cat("A: ",pairInfo$anode, " B: ",pairInfo$bnode, " C: ", pairInfo$cnode,"\n")
        cat("Pattern: ",pairInfo$adeco, " - ", pairInfo$bdeco1, " - ", 
            pairInfo$bdeco2, " - ", pairInfo$cdeco,"\n")
      }
      
      # gets involved nodes
      involvedNodes <- list(pairInfo$anode,pairInfo$bnode,pairInfo$cnode)
      inducedEdges <- private$getInducedSubgraph(involvedNodes)
      
      # checks the presence of an edge between anode and cnode
      acnode <- private$checkEdgeInArgumentEdges(pairInfo$anode, pairInfo$cnode, inducedEdges)
      
      # keeps on testing if there is no such edge
      if (acnode == FALSE){
        # now check if B belongs to AC separator
        # b should not belongs to ac separators
        acSeparator <- private$separators[[paste(pairInfo$anode,pairInfo$cnode,sep="")]]
        
        if (private$debug){
          cat("Separator set between ",pairInfo$anode," and ", pairInfo$cnode,"\n")
          print(acSeparator)
        }
        
        # If everything is ok make the change in the edges
        if (pairInfo$bnode %in% acSeparator){
          flag <- TRUE
        }
      } 
      
      if (private$debug){
        cat("\n -------------------------- checkConditionForStep14 end ---------------------\n")
      }
      
      # return flag
      return(flag)
    },
    
    #'##############################################################
    #' check if a pair of edges has the conditions for applying step 15
    #' arguments:
    #' @param pairInfo info about the pair of edges to consider
    #'##############################################################
    checkConditionForStep15 = function(pairInfo){
      flag <- FALSE
      
      if (private$debug){
        cat("\n -------------------------- checkConditionForStep15 begin ---------------------\n")
        cat("A: ",pairInfo$anode, " B: ",pairInfo$bnode, " C: ", pairInfo$cnode,"\n")
        cat("Pattern: ",pairInfo$adeco, " - ", pairInfo$bdeco1, " - ", 
            pairInfo$bdeco2, " - ", pairInfo$cdeco,"\n")
      }
            
      # now check if there is a wildcard-wildcard edge between A nd C
      egdes <- private$edgesInfo[pairInfo$anode == pairInfo$cnode,] 
      
      if (private$debug){
        cat("Edges between ",pairInfo$anode," and ", pairInfo$cnode," : ",nrow(edges),"\n")
      }
      
      # If everything is ok make the change in the edges
      if (nrow(edges) != 0){
        flag <- TRUE
        
        # change edge info: first edge remains and the second becomes to
        # block-wildcard
        edge <- changeEdgeInfo(pairInfo, 1, edgesInfo, "none", "none")
      }
      
      if (private$debug){
        cat("\n -------------------------- checkConditionForStep15 end ---------------------\n")
      }
      
      # return flag
      return(list(flag=flag, edge1=edge1, edge2=edge2))
    },
    
    #'##############################################################
    #' remove the edges in the path passed as argument. The edges
    #' are passed as argument in order to work in a recursive mode
    #' arguments:
    #' @param path path to consider
    #' @param edgesInfo data structure with edges information. This information
    #'                  will change as long as the method goes on
    #'##############################################################
    removeEdgesInPath = function(path, edgesInfo){
      # removes the edges of the path
      for(i in 1:length(path)){
        index1 <- i
        index2 <- i+1
        if (index2 > length(path)){
          index2=1
        }
        nodeA <- path[[index1]]
        nodeB <- path[[index2]]
        
        # gets the edge
        edge <- private$getEdgeInArgumentForNodes(nodeA, nodeB, edgesInfo)
        
        # removes the edge
        edgesInfo <- removeEdge(edge,edgesInfo)
      }
      
      # return edgesInfo
      return(edgesInfo)
    },
    
    #'##############################################################
    #' looks for an edge between nodeA and nodeB in the set of edges
    #' passed as argument
    #' arguments:
    #' @param nodeA first node of interest
    #' @param nodeB second node of interest
    #' @param edgersInfo data structure with edges information to analyze
    #'##############################################################
    getEdgeInArgumentForNodes = function(nodeA, nodeB, edgesInfo){
      if (private$debug){
        cat("\n -------------------------- getEdgeInArgumentoForNodes begin ---------------------\n")
        cat("NodeA: ",nodeA, "  nodeB: ",nodeB," \n")
        cat("List of edges: \n")
        print(edgesInfo)
      }
      
      # get edges having nodeA at one side and nodeB in the other
      edge <- edgesInfo[(edgesInfo$from == nodeA & edgesInfo$to == nodeB) |
                          (edgesInfo$to == nodeA & edgesInfo$from == nodeB), ]
      
      
      if (private$debug){
        cat("\n -------------------------- getEdgeInArgumentForNodes end ---------------------\n")
        cat("Edge: \n")
        print(edge)
        cat("-----------------------------------------------------------------------\n")
      }
      
      # make a common list with them
      return(edge)
    },
    
    #'##############################################################
    #' gets the edges related to the list of nodes passed as argument
    #' arguments:
    #' @param involvedNodes nodes to use for getting the induced subgraph
    #'##############################################################
    getInducedSubgraph = function(involvedNodes){
      # removes nodes with from or to not contained in involvedNodes
      involvedEdges <- private$edgesInfo[(private$edgesInfo$from %in% involvedNodes) &
                                           (private$edgesInfo$to %in% involvedNodes) ,]
      
      # return involvedEdges
      return(involvedEdges)
    },
    
    #'##############################################################
    #' checks if there is an edge between the nodes passed as first
    #' and second arguments and with the given decorations
    #' arguments:
    #' @param anode first node of interest
    #' @param bnode second node of interest
    #' @param aDeco decoration for first node
    #' @param bDeco decoration for second node
    #'##############################################################
    checkEdgeWithDecorations = function(anode, bnode, aDeco, bDeco){
      if (private$debug){
        cat("\n -------------------------- checkEdgeWithDecoration begin ---------------------\n")
        cat("Looking for direct edge between ",anode," and" ,bnode,"\n")
        cat("decoration for anode: ",aDeco," decoraction for bnode: ", bDeco,"\n")
      }
      
      flag <- FALSE
      
      # checks if there is such an edge
      edges <- private$edgesInfo[(private$edgesInfo$from == anode & private$edgesInfo$to == bnode) |
                             (private$edgesInfo$from == bnode & private$edgesInfo$to == anode),]
      
      if (private$debug){
        cat("Edges between the nodes: \n")
        print(edges)
      }
      
      # selects the edge and checks the decorations
      if (nrow(edges) != 0){
        edge <- edges[1,]
        
        # checks aDeco and bDeco
        if (edge$from == anode){
          leftDeco <- edge$left
          rightDeco <- edge$right
        }
        else{
          leftDeco <- edge$right
          rightDeco <- edge$left
        }
        
        if (private$debug){
          cat("Decoration for ",anode, " - ", leftDeco,"\n")
          cat("Decoration for ",bnode, " - ", rightDeco,"\n")
        }
        
        if (aDeco == "noblock"){
          flag <- (leftDeco != "block")
        }
        else{
          if (aDeco != "wildcard"){
            flag <- (aDeco == leftDeco)
          }
        }
        
        # the same for bdeco
        if (bDeco == "noblock"){
          flag <- (rightDeco != "block")
        }
        else{
          if (bDeco != "wildcard"){
            flag <- (bDeco == rightDeco)
          }
        }
      }
      
      if (private$debug){
        cat("Result of check: ",flag,"\n")
        cat("\n -------------------------- checkEdgeWithDecorations end ---------------------\n")
      }
      
      # returns flag
      return(flag)
    },
    
    #'##############################################################
    #' checks if there is an edge between the nodes passed as first
    #' and second arguments
    #' arguments:
    #' @param anode first node
    #' @param bnode second node
    #' @param inducedEdges info about edges to analyze
    #'##############################################################
    checkEdgeInArgumentEdges = function(anode, cnode, inducedEdges){
      if (private$debug){
        cat("\n -------------------------- checkEdgeInArgumentEdges begin ---------------------\n")
        cat("Looking for direct edge between ",anode," and" ,cnode,"\n")
        cat("Edges in induced graph\n")
        print(inducedEdges)
      }
      
      flag <- FALSE
      
      # checks if there is such an edge
      edge <- inducedEdges[(inducedEdges$from == anode & inducedEdges$to == cnode) |
                             (inducedEdges$from == cnode & inducedEdges$to == anode),]
      
      if (nrow(edge) != 0){
        flag <- TRUE
      }
      
      if (private$debug){
        cat("Result of check: ",flag,"\n")
        cat("\n -------------------------- checkEdgeInArgumentEdges end ---------------------\n")
      }
      
      # returns flag
      return(flag)
    },
    
    #'##############################################################
    #' changes the info of an edge with the information passed as first argument
    #' the second argument shows the edge to change: 1 or 2
    #' the information about edges is required in order to select the proper
    #' edges and making the change on them
    #' arguments: 
    #' @param  pairInfo info about the pair of edges containing the edge to change
    #' @param  index 1 or 2, showing the edge in pairInfo to change
    #' @param  edgesInfo information about edges. This info will be modified after
    #              the change
    #' @param  leftDeco new decoration for an extreme
    #' @param  rightDeco new decoration for another extreme
    #'##############################################################
    changeEdgeInfo = function(pairInfo, index, edgesInfo, leftDeco, rightDeco){
      if (private$debug){
        cat("\n -------------------------- changeEdgeInfo begin ---------------------\n")
        cat("Info about pair: \n")
        print(pairInfo)
        cat("Edge to change: ",index,"\n")
        cat("New decoration to apply: ",leftDeco," - ", rightDeco, "\n")
      }
      
      # selects the nodes of the edge to change
      reverted <- FALSE
      if (index == 1){
        if (pairInfo$e1reverted == FALSE){
          from <- pairInfo$anode
          to <- pairInfo$bnode
        }
        else{
          from <- pairInfo$bnode
          to <- pairInfo$anode
          reverted <- TRUE
        }
      }
      else{
        # information about edge number 2
        if (pairInfo$e2reverted == FALSE){
          from <- pairInfo$bnode
          to <- pairInfo$cnode
        }
        else{
          from <- pairInfo$cnode
          to <- pairInfo$bnode
          reverted <- TRUE
        }
      }
      
      if (private$debug){
        cat("In original information: from = ",from, " - to = ",to," reverted: ",reverted, "\n")
      }
      
      # now gets the edge from edgesInfo
      edge <- edgesInfo[(edgesInfo$from == from & edgesInfo$to == to),]
      
      if (private$debug){
        cat("Original edge: \n")
        print(edge)
      }
      
      # now change the decorations according to the reversion
      if (reverted){
        temp <- leftDeco
        leftDeco <- rightDeco
        rightDeco <- temp
      }
      
      # change the decorations
      if(leftDeco != "wildcard"){
        edge$left=leftDeco
      }
      if (rightDeco != "wildcard"){
        edge$right=rightDeco
      }
      
      if (private$debug){
        cat("Final edge: \n")
        print(edge)
        cat("\n -------------------------- changeEdgeInfo end ---------------------\n")
      }
      
      # return the edge
      return(edge)
    },
    
    #'##############################################################
    #' method for changing the decorations of an edge between anode
    #' and bnode
    #' arguments:
    #' @param anode first node
    #' @param bnode second node
    #' @param anodeDeco deco for first node extreme
    #' @param bnodedeco deco for second node extreme
    #'##############################################################
    changeEdgeDecoration = function(anode, bnode, anodeDeco, bnodeDeco){
      if (private$debug){
        cat("\n -------------------------- changeEdgeDecoration begin ---------------------\n")
        cat("Changing edge between ",anode," and ",bnode, " to decos ",anodeDeco, " - ", bnodeDeco,"\n")
      }
      
      changed <- FALSE
      edge <- NULL
      # looks for anode - bnode edge (with from - to direction)
      edge <- private$edgesInfo[private$edgesInfo$from == anode &
                                private$edgesInfo$to == bnode,]
      index <- which(private$edgesInfo$from == anode & private$edgesInfo$to == bnode)
      
      # if this edge exists, just change the decoration
      if(length(index) == 1){
        if (anodeDeco != "wildcard"){
          edge$left <- anodeDeco
        }
        if (bnodeDeco != "wildcard"){
          edge$right <- bnodeDeco
        }        
        
        # if working in pc mode avoid block-block edges
        blockBlock <- FALSE
        if (edge$left == "block" & edge$right == "block"){
          blockBlock <- TRUE
        }
        
        # change the edge if needed
        if (!private$pc | (private$pc & !blockBlock)){
          private$edgesInfo[index,] <- edge
          changed <- TRUE
        }
      }
      else{
        # the edge must be in to - from orientation
        edge <- private$edgesInfo[private$edgesInfo$to == anode &
                                    private$edgesInfo$from == bnode,]
        index <- which(private$edgesInfo$to == anode & private$edgesInfo$from == bnode)
        
        # change the edge decoration if it exists
        if(length(index) == 1){
          if (bnodeDeco != "wildcard"){
            edge$left <- bnodeDeco
          }
          if(anodeDeco != "wildcard"){
            edge$right <- anodeDeco
          }
          # if working in pc mode avoid block-block edges
          blockBlock <- FALSE
          if (edge$left == "block" & edge$right == "block"){
            blockBlock <- TRUE
          }
          
          # change the edge if needed
          if (!private$pc | (private$pc & !blockBlock)){
            private$edgesInfo[index,] <- edge
            changed <- TRUE
          }
        }
      }
      if (private$debug){
        cat("changed: ",changed,"\n")
        cat("\n -------------------------- changeEdgeDecoration end ---------------------\n")
      }
      
      # return changed
      return(changed)
    },
    
    #'##############################################################
    #' method for checking the existance of a loop between origin and dest
    #' without chords and using edges with a certain decoration
    #' arguments:
    #' @param origin first node to consider in the path
    #' @param dest last node in the path
    #' @param originDeco decoration for origin side in the edge starting the path
    #' @param destDeco decoration for dest in the edge ending the path
    #' @param path path to analyze
    #' @param edges info about edges in the model
    #'##############################################################
    isReachableWithLoopAndDeco = function(origin, dest, originDeco, destDeco, path, edges){
      if (private$debug){
        cat("\n -------------------------- isReachableWithLoopAndDeco begin ---------------------\n")
        cat("Origin: ",origin, " destination: ", dest, " path: ")
        cat("Origin deco: ", originDeco, " dest deco: ", destDeco)
        print(path)
        cat("Edges to consider: ")
        print(edges)
      }
      
      # initializes an empty result
      result <- list("path"=path,"flag"=FALSE,"length"=length(path),"loop"=FALSE)
      
      # base case: dest reached and non empty path
      if (origin == dest & length(path) != 0){
        if (private$debug){
          cat("Base case: dest and origin are the same. Path length: ",length(path),"\n")
        }
        
        # sets flag to true
        result$flag <- TRUE
        result$loop <- TRUE
        
        # in any case return
        return(result)
      }
      
      # if this point is reached, it is needed to keep on searching
      # get adjacents of current node  but taking into account edges 
      # without blocks
      candidateEdges <- private$getEdgesInArgumentForNode(origin, originDeco, destDeco, edges)
      
      if (private$debug){
        cat("Obtained the list of adjacent to ",origin, "\n")
        cat("Candidate edges to consider: ",nrow(candidateEdges),"\n")
      }
      
      # considers every candidate
      if (private$containsData(candidateEdges)){
        for(i in 1:nrow(candidateEdges)){
          # selects the candidate
          candidate <- candidateEdges[i,]
          
          if (private$debug){
            cat("Candidate edge under analysis: \n")
            print(candidate)
          }
          
          # gets the following node
          other <- private$getFollowingNode(origin, candidate)
          
          if (private$debug){
            cat("New destination: ",other,"\n")
          }
          
          # if candidate is in path there is no need to go on. In any other
          # case go on analyzing edges
          if (!(other %in% path)){
            if (private$debug){
              cat("Analyzing edge ", candidate$from, " - ", candidate$to, "\n")
            }
            
            # adds other to path
            path <- append(path,other)
            
            # removes the edge
            edgesNew <- private$removeEdgeInArgument(candidate, edges)
            
            if (private$debug){
              cat("Actual path from origin: \n")
              print(path)
            }
            
            # and checks if the other dest (the final destination, passed
            # as argument) can be reached from destNode
            result <- private$isReachableWithLoopAndDeco(other, dest , originDeco, destDeco, 
                                                 path, edgesNew)
            
            # if the result is true, no need to keep on searching
            if(result$flag == TRUE){
              break
            }
            else{
              # removes the node under consideration from path
              path <- path[path != other]
              if (private$debug){
                cat("Removing ",other," from path\n")
                print(path)
              }
            }
          }
          else{
            # the node to consider already belongs to path
            result$flag <- FALSE
            result$loop <- TRUE
            result$path <- append(path,other)
            break;
          }
        }
        
        if (private$debug){
          cat("\n -------------------------- isReachableWithLoopAndDeco end ---------------------\n")
          cat("flag: ",result$flag,"\n")
          cat("loop: ",result$loop,"\n")
          cat("-------------------------------------------------------------------\n")
        }
      }
      
      # return result
      return(result)
    },
    
    #'##############################################################
    #' gets the node at the opposite extreme of origin in edge
    #' arguments:
    #' @param origin origin node
    #' @param edge edge of interest
    #'##############################################################
    getFollowingNode = function(origin, edge){
      if (origin == edge$to){
        other <- edge$from
      }
      else{
        other <- edge$to
      }
      
      # return other
      return(other)
    },
    
    #'##############################################################
    #' removes the information related to a certain edge but not on the
    #' data member but on the argument
    #' arguments:
    #' @param edge edge to remove
    #' @param edgesInfo list of edges where the edge will be removed
    #'##############################################################
    removeEdgeInArgument = function(edge, edgesInfo){
      # removes the used link from edgesInfo
      ind <- which(with(edgesInfo, edgesInfo$from == edge$from & 
                          edgesInfo$to == edge$to))
      edgesInfo <- edgesInfo[ -ind, ]
      return(edgesInfo)
    },
    
    #'##############################################################
    #' method for getting quartets of nodes in order to apply rule 4
    #' arguments:
    #' @param pairs pairs of candidate pairs of edges
    #'##############################################################
    getQuartets = function(pairs){
      if (private$debug){
        cat("\n -------------------------- getQuartets begin ---------------------\n")
        cat("Number of pairs: ",nrow(pairs),"\n")
        cat("Pairs to consider: \n")
        cat(class(pairs),"\n")
        print(pairs)
        cat("\n")
      }
      candidatePairs <- NULL
      
      # considers every pair in restpairs
      index <- seq(1:(nrow(pairs)-1))
      
      # form candidate pairs of edges: only if there is more than
      # a edge
      if (nrow(pairs) > 1){
        candidatePairs <- lapply(index, private$getQuartetsForPair, df=pairs)
        candidatePairs <- private$delete.NULLs(candidatePairs)
      }
      
      if (private$debug){
        cat("\n -------------------------- getQuartets end ---------------------\n")
      }
      return(candidatePairs)
    },
    
    #'##############################################################
    #' get quartets for a given pair of edges
    #' arguments:
    #' @param index edge considered as base 
    #' @param df data estructure with information to consider
    #'##############################################################
    getQuartetsForPair = function(index, df){
      # selects the base pair
      basepair <- df[index,]
      
      # gets the rest of elements in pairs
      restpairs <- df[(index+1):nrow(df),]
      
      # considers every pair in restpairs
      index <- seq(1:nrow(restpairs))
      
      # now checks if basepair and some of the pairs in restpairs
      # can be considered for applying rule 4. It is needed a and
      # c nodes must match between them
      quartets <- lapply(index,function(index){
        secondpair <- restpairs[index,]
        if ((secondpair$anode == basepair$anode) & (secondpair$cnode == basepair$cnode)){
          return(list("pair1"=basepair,"pair2"=secondpair))
        }
      })
      
      # removes nulls
      quartets[sapply(quartets, is.null)] <- NULL
      return(quartets)
    },
    
    #'##############################################################
    #' remove the edges in the path passed as argument
    #' arguments:
    #' @param path path to consider
    #' @param edgesInfo data structure with edges info
    #'##############################################################
    removeEdgesInArgumentInPath = function(path, edgesInfo){
      if (private$debug){
        cat("\n -------------------------- removeEdgesInArgumentInPath begin ---------------------\n")
        cat("Path: \n")
        print(path)
      }
      # removes the edges of the path
      for(i in 1:length(path)){
        index1 <- i
        index2 <- i+1
        if (index2 > length(path)){
          index2=1
        }
        nodeA <- path[[index1]]
        nodeB <- path[[index2]]
        
        # gets the edge
        edge <- private$getEdgeInArgumentForNodes(nodeA, nodeB, edgesInfo)
        
        # removes the edge
        edgesInfo <- private$removeEdgeInArgument(edge,edgesInfo)
      }
      
      if (private$debug){
        cat("\n -------------------------- removeEdgesInArgumentInPath end ---------------------\n")
      }
      
      # return edgesInfo
      return(edgesInfo)
    },
    
    #'##############################################################
    #' detects duplicated nodes in a path and remove the
    #' links related to the cycle
    #' arguments: 
    #' @param path: path to consider
    #' @param edges: set of edges to analyze
    #'##############################################################
    removeEdgesInArgumentInCycle = function(path, edges){
      if (private$debug){
        cat("\n -------------------------- removeEdgesInArgumentInCycle begin ---------------------\n")
        cat("Path with cycle: \n")
        print(path)
      }
      
      # gets the duplicated node
      cycleNode <- path[duplicated(path)]
      
      if (private$debug){
        cat("Node in cycle origin: ",cycleNode[[1]])
      }
      
      # gets the first index containing 
      repetitions <- which(with(path,path == cycleNode[[1]]))
      index1 <- repetitions[[1]]
      index2 <- repetitions[[2]]
      
      if (private$debug){
        cat("Cycle between positions ",index1," and ",index2)
      }
      
      # gets the path from index1 to index2
      cycle <- path[index1:index2]
      
      # now remove the edges in the cycle
      edges <- private$removeEdgesInArgumentInPath(cycle,edges)
      
      if (private$debug){
        cat("\n -------------------------- removeEdgesInArgumentCycle begin ---------------------\n")
      }
      
      # return edges
      return(edges)
    },
    
    #'##############################################################
    #' method to filter the pair of edges according to the separator
    #' set of the extreme nodes. The boolean flag shows if the common
    #' node must be contained on the the separator set or not
    #' arguments:
    #' @param pairs: pairs of edges
    #' @param commonIncluded: boolean flag to state if the common node
    #'                   must be contained in separators or not
    #'##############################################################
    filterPairsBySeparators = function(pairs, commonIncluded){
      
      # considers every pair in pairs
      index <- seq(1:(nrow(pairs)-1))
      
      # form candidate pairs of edges: only if there is more than
      # a edge
      if (nrow(pairs) > 1){
        pairs <- lapply(index, private$checkPairSeparators, df=pairs, commonIncluded=commonIncluded)
        # remove nulls
        pairs <- private$delete.NULLs(pairs)
      }
           
      result <- NULL

      if (private$containsData(pairs)){
        result <- data.table::rbindlist(pairs) 
      }

      # return the filtered pairs
      return(result)
    },
    
    #'##############################################################
    #' checks if the node in common belongs or not to the
    #' separator set
    #' arguments:
    #' @param index: index of pairs to consider
    #' @param df: data structure with pairs of edges
    #' @param commonIncluded: boolean flag
    #'##############################################################
    checkPairSeparators = function(index, df, commonIncluded){
      pair <- df[index,]
      # get the nodes in both edges
      nodes1 <- c(pair$from.x, pair$to.x)
      nodes2 <- c(pair$from.y, pair$to.y)
      
      # gets B (common node)
      bNode <- intersect(nodes1, nodes2)
     
      # gets the nodes at the extremes
      diff <- setdiff(union(nodes1,nodes2),bNode)
      aNode <- diff[1]
      cNode <- diff[2]
      
      # gets the separators for aNode and cNode
      acSeparator <- private$separators[[paste(aNode,cNode,sep="")]]
      
      # initialized FLAG to false and checks if b belongs to acSeparator
      contained <- (bNode %in% acSeparator)
      outputPair <- NULL
      if ((commonIncluded & contained)  | (!commonIncluded & !contained)){
        outputPair <- pair
      }
    
      # return outputPair
      return(outputPair)
    },
    
    #'##############################################################
    #' checks if there is an edge for two nodes with a certain pattern of
    #' decorations:
    #' arguments: 
    #' @param anode first node of interest
    #' @param bnode second node of interest
    #' @param leftDeco decoration for anode
    #' @param rightDeco decoration for bnode
    #'##############################################################
    checkEdgeForNodesWithPattern = function(anode, bnode, leftDeco, rightDeco){
      if (private$debug){
        cat("\n -------------------------- checkEdgeForNodesWithPattern begin ---------------------\n")
      }
      directedEdges <- private$getEdgesWithPatternFromTo(leftDeco, rightDeco)
      if (!is.null(directedEdges)){
        directedEdges <- directedEdges[directedEdges$from == anode & directedEdges$to == bnode,]
      }
      
      invertedEdges <- private$getEdgesWithPatternFromTo(leftDeco,rightDeco)
      if (!is.null(invertedEdges)){
        invertedEdges <- invertedEdges[invertedEdges$from== bnode & invertedEdges$to == anode,]
      }
      
      flag <- ((!is.null(directedEdges) & nrow(directedEdges) != 0) |
               (!is.null(invertedEdges) & nrow(invertedEdges) != 0))
        
        if (private$debug){
          cat("\n -------------------------- checkEdgeForNodesWithPattern end ---------------------\n")
        }
      
      return(flag)
    },
    
    #'##############################################################
    #' get the set of edges with a certain pattern in decorations
    #' This method looks into the set of edges defined in the data
    #' member edgesInfo
    #' arguments: 
    #' @param fromDeco decoration for from side
    #' @param toDeco  decoration for to side
    #'##############################################################
    getEdgesWithPatternFromTo = function(fromDeco, toDeco){
      if (private$debug){
        cat("\n -------------------------- getEdgesWithPatternFromTo begin ---------------------\n")
        cat("fromDeco: ",fromDeco," toDeco: ",toDeco,"\n")
      }
      
      # select edges with the goven decoration on the left side
      if (fromDeco == "noblock"){
        edges <- private$edgesInfo[private$edgesInfo$left != "block",]
      }
      else{
        if (fromDeco != "wildcard"){
          edges <- private$edgesInfo[private$edgesInfo$left == fromDeco,]
        }
      }
      
      # the same for the right decoration
      if (toDeco == "noblock"){
        edges <- edges[edges$right != "block",]
      }
      else{
        if (toDeco != "wildcard"){
          edges <- edges[edges$right == toDeco,]
        }
      }
      
      if (private$debug){
        cat("\n -------------------------- getEdgesWithPatternFromTo end ---------------------\n")
        cat("Returned edges: \n")
        print(edges)
      }
      
      # return edges
      return(edges)
    },
    
    #'##############################################################
    #' updates the graph according to the set of edges
    #' With his method graph data member is updated to match the
    #' information contained into edgesInfo data member
    #'##############################################################
    updateGraph = function(){
      
      # considers every edge in edgesInfo
      for(i in 1:nrow(private$edgesInfo)){
        edge <- private$edgesInfo[i,]
        if (edge$left == "arrow" & edge$right == "none"){
          # sets the orientation in the graph
          private$graph <- bnlearn::set.arc(private$graph,from=edge$to,to=edge$from,check.cycles=FALSE,debug=FALSE)
        }
        else{
          if (edge$right == "arrow" & edge$left == "none"){
            # sets the orientation in the graph
            private$graph <- bnlearn::set.arc(private$graph,from=edge$from,to=edge$to,check.cycles=FALSE,debug=FALSE)
          }
        }
      }
    },
    
    #'##############################################################
    #' gets vstructures with mampc perspective
    #' arguments:
    #' @param edges list of edges to analyze
    #'##############################################################
    getModelVStructuresDataFrame = function(edges){
      if (private$debug){
        cat("\n -------------------------- getModelVStructuresDataFrame begin ---------------------\n")
      }
      
      # creates the data frame structure for storing v-structures
      vsdf <- data.frame(X=character(), Z=character(), 
                         Y=character(), type=numeric(), stringsAsFactors=FALSE)
      
      # get pairs with a common node as first step
      pairs <- private$matchEdgesWithCommonNode(edges)

      # # generate a data frame if there are pairs; in any other case
      # NULL will be returned
      candidatePairsdf <- private$generateAlignmentDataFrame(pairs)

      # remove candidates according to moral parameter
      if (private$moral == FALSE){
        toRemove <- c()
        if (!is.null(candidatePairsdf) && nrow(candidatePairsdf) != 0){
          for(i in 1:nrow(candidatePairsdf)){
            # selects the pair
            pair <- candidatePairsdf[i,]

            # check the condition
            if (pair$anode %in% private$getAdjacents(pair$cnode) | 
                pair$cnode %in% private$getAdjacents(pair$anode)){
              toRemove <- c(toRemove, i)
            }
          }
        }

        # remove pairs where anode and cnode are adjacents
        if (length(toRemove) != 0){
          candidatePairsdf <- candidatePairsdf[-toRemove,]
        }
      }
      
      # checks if there a pair and selects pairs of edges with the required pattern
      if (private$containsData(candidatePairsdf)){
          
        # selects the pattern block - wildcard - wildcard - wildcard
        vs <- private$selectPairPattern(candidatePairsdf, "none","arrow", 
                                          "arrow", "none")
          
        # inserts the information into vsdf
        vsdf <- private$insertVSInDataFrame(vsdf, vs, 1)
          
        # selects the pattern none - arrow - none -none
        vs <- private$selectPairPattern(candidatePairsdf, "none","arrow", 
                                          "none", "none")
          
        # inserts the information into vsdf
        vsdf <- private$insertVSInDataFrame(vsdf, vs, 2)
          
        # the same for patter none - arrow - arrow - arrow
        vs <- private$selectPairPattern(candidatePairsdf, "none","arrow", 
                                          "arrow", "arrow")
        
        # inserts the information into vsdf
        vsdf <- private$insertVSInDataFrame(vsdf, vs, 3)
          
        # patter none - none - arrow - none
        vs <- private$selectPairPattern(candidatePairsdf, "none","none", 
                                          "arrow", "none")
          
        # inserts the information into vsdf
        vsdf <- private$insertVSInDataFrame(vsdf, vs, 4)
          
        # patter none - none - arrow - arrow
        vs <- private$selectPairPattern(candidatePairsdf, "none","none", 
                                          "arrow", "arrow")
          
        # inserts the information into vsdf
        vsdf <- private$insertVSInDataFrame(vsdf, vs, 5)
          
        # patter arrow - arrow - none -none
        vs <- private$selectPairPattern(candidatePairsdf, "arrow","arrow", 
                                          "none", "none")
          
        # inserts the information into vsdf
        vsdf <- private$insertVSInDataFrame(vsdf, vs, 6)
          
        # pattern arrow - arrow - arrow - none
        vs <- private$selectPairPattern(candidatePairsdf, "arrow","arrow", 
                                          "arrow", "none")
          
        # inserts the information into vsdf
        vsdf <- private$insertVSInDataFrame(vsdf, vs, 7)
          
        # pattern arrow - arrow - arrow - arrow
        vs <- private$selectPairPattern(candidatePairsdf, "arrow","arrow", 
                                          "arrow", "arrow")
          
        # inserts the information into vsdf
        vsdf <- private$insertVSInDataFrame(vsdf, vs, 8)
      }
      
      if (private$debug){
        cat("\n -------------------------- getModelVStructuresDataFrame end ---------------------\n")
      }
      
      # return vs
      return(vsdf)
    },
    
    #'##############################################################
    #' function to process a pair of edges checking if they compose
    #' a v-structure
    #' arguments:
    #' @param pair: pair of edges to consider (as indices)
    #' @param edges: list of edges to analyze
    #'##############################################################
    checkVStructure = function(pair, edges){
      edge1 <- edges[pair[1],]
      edge2 <- edges[pair[2],]
      vstructure <- NULL
      if (edge1["to"] == edge2["to"]){
        vstructure <- list("A"=edge1[["from"]], "B"=edge1[["to"]], "C"=edge2[["from"]])
      }
      return(vstructure)
    },
    
    #'##############################################################
    #' insert vstructures stored in pairs in the data frame passed
    #' as a first argument, avoiding repetitions
    #' arguments:
    #' @param vsdf: data estructure to return
    #' @param pairs: pair of edges to consider
    #' @param type: argument to be passes to the auxiliar method 
    #'##############################################################
    insertVSInDataFrame = function(vsdf, pairs, type){      
      # checks if there is something to insert
      if (private$containsData(pairs)){
        for(i in 1:nrow(pairs)){
          pair <- pairs[i,]
          vsdf <- private$checkAndInsertVS(pair, vsdf, type)
        }
      }
      
      # return vsdf
      return(vsdf)
    },
    
    #'##############################################################
    #' check is a vstructure is contained in the data frame
    #' arguments:
    #' @param pair pair of edges to consider
    #' @param vsdf data structure where information must be stored
    #' @param type type of v-structure
    #'##############################################################
    checkAndInsertVS = function(pair, vsdf, type){
      insertPosition <- nrow(vsdf)+1
      
      # composes the natural vstructure info
      vs <- list(X=pair$anode, Z=pair$bnode, Y=pair$cnode, type=type)
      vsinverted <- list(X=pair$cnode, Z=pair$bnode, Y=pair$anode, type=type)
      # checks if present
      if (!private$find.row(vsdf, vs) & !private$find.row(vsdf,vsinverted)){
        vsdf[insertPosition,] <- vs
      }
      
      # return df
      return(vsdf)
    },
      
    #'##############################################################
    #' find duplicates of a row in a data frame
    #' arguments:
    #' @param df: data frame to consider
    #' @param row: row to consider
    #'##############################################################
    find.row = function(df, row){
      if (private$debug){
        cat("\n -------------------------- find.row begin ---------------------\n")
        cat("row to look for: \n")
        print(row)
        cat("------------------------------------------------------------------------------")
      }
      
      # initializes found to FALSE
      found <- FALSE
      
      if (nrow(df) != 0){
        matched <- df[(df$X == row$X & df$Y == row$Y & df$Z == row$Z), ]
        found <- (nrow(matched) != 0)
      }
      
      if (private$debug){
        cat("\n -------------------------- find.row end ---------------------\n")
        cat("Returning ", found, "\n")
        cat("------------------------------------------------------------------------------")
      }
      
      return(found)
    },
    
    #'##############################################################
    #' compares vstructures passed as argument
    #' arguments: 
    #' @param targetVsdf: info about true v-structures
    #' @param currentVsdf: info about learnt v-structures
    #'##############################################################
    compareVStructures = function(targetVsdf, currentVsdf){
      if (private$debug){
        cat("\n -------------------------- compareVStructures begin ---------------------\n")
        cat("Edges in current structure: \n")
        print(private$edgesInfo)
        cat("V-structures in target structure: \n")
        print(targetVsdf)
        cat("V-structures in current structure\n")
        print(currentVsdf)
        cat("------------------------------------------------------------------------------")
      }
      # initializes variables
      # tp: present in both sets
      # fp: present in current and missing in target
      # fn: present in target and missing in current
      tp <- 0
      fp <- 0
      fn <- 0
      
      # considers structures in targetVsdf and check if present in currentVsdf
      if (private$debug){
        cat("Comparing target against current: \n")
      }
      
      if(nrow(targetVsdf) != 0){
        for(i in 1:nrow(targetVsdf)){
          vs <- targetVsdf[i,]
          vsInverted <- list(X=vs$Y,Z=vs$Z,Y=vs$X,type=vs$type)
        
          if (private$debug){
            cat("Checking : \n")
            print(vs)
            cat(" and \n")
            print(vsInverted)
          }
          if (private$find.row(currentVsdf,vs)){
            tp <- tp+1
            if (private$debug){
              cat("Match in normal direction: ",tp,"\n")
            }
          }
          else{
            # checks the opposite direction
            if (private$find.row(currentVsdf,vsInverted)){
              tp <- tp+1
              if (private$debug){
                cat("Match in opposite direction: ",tp,"\n")
              }
            }
            else{
              fn <- fn+1
              if (private$debug){
                cat("Not-found v-structure: ",fn,"\n")
              }
            }
          }
        }
      }
      
      if (private$debug){
        cat("\nComparing current against target: \n")
      }
      
      # considers structures in currentVsdf and check if missing in targetVsdf
      if (nrow(currentVsdf) != 0){
        for(i in 1:nrow(currentVsdf)){
          vs <- currentVsdf[i,]
          vsInverted <- list(X=vs$Y,Z=vs$Z,Y=vs$X,type=vs$type)
          if (!private$find.row(targetVsdf,vs) & !private$find.row(targetVsdf,vsInverted)){
            fp <- fp+1
            if (private$debug){
              cat("\nV-structure not found on target structure: ",fn,"\n")
            }
          }
        }
      }
      
      # now computed the corresponding measures
      if (private$debug){
         cat("COUNTERS: tp: ",tp," fn: ", fn, " fp: ",fp,"\n")
      }
      recall <- tp/(tp+fn)
      precision <- tp/(tp+fp)
      
      if (private$debug){
        cat("Recall: ",recall,"\n");
        cat("Precision: ", precision,"\n")
        cat("\n -------------------------- compareVStructures end ---------------------\n")
      }
      
      # return the results
      return(list(recall=recall,precision=precision))
    },
    
    #'##############################################################
    #' delete null values from a list
    #' arguments:
    #' @param x.list: list to consider
    #'##############################################################
    delete.NULLs  =  function(x.list){   # delele null/empty entries in a list
      x.list[unlist(lapply(x.list, length) != 0)]
    },

    #'##############################################################
    #' checks the data is not null and contains data
    #' arguments:
    #' @param data: data structure to consider
    #'##############################################################
    containsData = function(data){
      flag <- FALSE
      
      if (!is.null(data)){
        flag <- TRUE
        if((class(data)[1] == "list")){
          flag <- (length(data) != 0)
        }
        else{
          flag <- (nrow(data) != 0)
        }
      }
      return(flag)
    }
  )
)

#'##############################################################
#' easy constructor of MampcSearch class. 
#' NOTE: Observe default values for some parameters
#'##############################################################
buildObject <- function(net, data, moral=TRUE, pc=FALSE, alpha=0.05, edges=NULL, debug=FALSE){
  MAMPCGSearch$new(net, data, moral, pc, alpha, edges, debug)
}

