
library(glue)

hammingDistance <- function(x1,x2) sum(x1 != x2)

#number of nodes in a network loaded with the load() function
numberOfGenes <- function(net){
    return (length(net[[1]])) # We read from the imported network its number of nodes
}

#ATTRACTORS TO FILE
saveAttractors <- function(attractors, prefixName, path="."){
    for (i in seq(1,length(attractors$attractors))){
        att <- getAttractorSequence(attractors, i)
        write.table(att, glue('{path}/{prefixName}_att_{i}.csv'), sep = ",", col.names = colnames(att),row.names=FALSE)
    }
}

fromIntegerToBitVector <- function(integer, numOfBits) {rev(as.integer(intToBits(integer))[1:numOfBits])}
#Generate a RANDOM Boolean state of a specified length
rndBooleanState     <- function(vectorLength){ return (sample(c(0,1), vectorLength, replace=TRUE)) }
rndBooleanVector    <- function(vectorLength, bias){ return (sample(c(0,1), vectorLength, replace=TRUE, prob = c(1-bias, bias))) }

# From Boolean vector to Logical expression with minterms
sumOfProducts <- function(gene, incomingGenes, booleanFunction){
    numGenes <- length(incomingGenes)
    string <- glue('Gene{gene}, ')
    atLeastOnePositiveOutput <- FALSE
    #print(booleanFunction)
    if (length(unique(booleanFunction)) > 1){ #caso in cui non siano tutti 0 o tutti 1
        for (index in seq(1,length(booleanFunction))){
            if (booleanFunction[index] == 1){
                atLeastOnePositiveOutput <- TRUE
                #print(index)
                #print(booleanFunction[index])
                minterm <- rev(as.numeric(intToBits(index-1))[1:numGenes])
                #print(minterm)
                string <- paste0(string, '(')
                for (i in seq(1,length(minterm))){
                    if (minterm[i] == 1){
                        string <- paste0(string, glue('Gene{incomingGenes[i]}'))
                    } else {
                        string <- paste0(string, glue('!Gene{incomingGenes[i]}'))
                    }
                    string <- paste0(string, ' & ')
                }
                string <- gsub('.{3}$', '', string)
                string <- paste0(string, ')')
                string <- paste0(string, ' | ')
            }
        }
        if (atLeastOnePositiveOutput){
            string <- gsub('.{3}$', '', string)
        }
    } else {
        string <- paste0(string, glue('{booleanFunction[1]}'))
    }
    return(string)
}

rbn <- function(n, k, p, selfLoops=FALSE){
    bn <- vector("list", length = n)
    nodes <- c(1:n)
    for (i in (seq(1,n))){
        bn[[i]][["id"]]    <- i
        bn[[i]][["k"]]          <- k
        bn[[i]][["p"]]          <- p
        if (selfLoops) {
            bn[[i]][["incoming"]]   <- sample(nodes, k, replace=FALSE)
        } else {
            bn[[i]][["incoming"]]   <- sample(nodes[!nodes == i], k, replace=FALSE)
        }
        bn[[i]][["func"]]       <- rndBooleanVector(2^bn[[i]][["k"]], p)
        bn[[i]][["expr"]]       <- sumOfProducts(i, bn[[i]][["incoming"]], bn[[i]][["func"]])
    }
    return(bn)
} 

#Draw an RBN from the subset defined by the setOfAllowedFunctions
rbnSubset <- function(n, k, p, selfLoops=FALSE, setOfAllowedFunctions){
    bn <- vector("list", length = n)
    nodes <- c(1:n)
    for (i in (seq(1,n))){
        bn[[i]][["id"]]    <- i
        bn[[i]][["k"]]          <- k
        bn[[i]][["p"]]          <- p
        if (selfLoops) {
            bn[[i]][["incoming"]]   <- sample(nodes, k, replace=FALSE)
        } else {
            bn[[i]][["incoming"]]   <- sample(nodes[!nodes == i], k, replace=FALSE)
        }
        bn[[i]][["func"]]       <- sample(setOfAllowedFunctions, 1, replace=TRUE)[[1]]
        bn[[i]][["expr"]]       <- sumOfProducts(i, bn[[i]][["incoming"]], bn[[i]][["func"]])
    }
    return(bn)
} 

saveNetworkToFile <- function(bn, filename){
    head <- "targets, factors"
    l <- lapply(bn, `[[`, "expr")
    l <- paste(unlist(l), collapse="\n")
    l <- paste(head, l, sep="\n")
    write(l,file=filename)
}

################
###### OR ######
################
constORSelfLoop <- function(bn, gene){
    previousIncoming            <- bn[[gene]][["incoming"]]
    previousIncoming[1]         <- gene
    bn[[gene]][["incoming"]]    <- previousIncoming
    bn[[gene]][["func"]]        <- c(0 ,rep(1, 2^bn[[gene]][["k"]] - 1))
    bn[[gene]][["expr"]]        <- sumOfProducts(gene, bn[[gene]][["incoming"]], bn[[gene]][["func"]])
    return(bn)
}
augmORSelfLoop <- function(bn, gene){
    previousIncoming            <- bn[[gene]][["incoming"]]
    previousIncoming            <- c(gene, previousIncoming)
    previousOutputLength        <- length(bn[[gene]][["func"]])
    bn[[gene]][["k"]]           <- bn[[gene]][["k"]] + 1
    bn[[gene]][["incoming"]]    <- previousIncoming
    bn[[gene]][["func"]]        <- c(bn[[gene]][["func"]], rep(1, previousOutputLength))
    bn[[gene]][["expr"]]        <- sumOfProducts(gene, bn[[gene]][["incoming"]], bn[[gene]][["func"]])
    return(bn)
}

################
###### AND #####
################
constANDSelfLoop <- function(bn, gene){
    previousIncoming            <- bn[[gene]][["incoming"]]
    previousIncoming[1]         <- gene
    bn[[gene]][["incoming"]]    <- previousIncoming
    bn[[gene]][["func"]]        <- c(rep(0, 2^bn[[gene]][["k"]] - 1), 1)
    bn[[gene]][["expr"]]        <- sumOfProducts(gene, bn[[gene]][["incoming"]], bn[[gene]][["func"]])
    return(bn)
}
augmANDSelfLoop <- function(bn, gene){
    previousIncoming            <- bn[[gene]][["incoming"]]
    previousIncoming            <- c(gene, previousIncoming)
    previousOutputLength        <- length(bn[[gene]][["func"]])
    bn[[gene]][["k"]]           <- bn[[gene]][["k"]] + 1
    bn[[gene]][["incoming"]]    <- previousIncoming
    bn[[gene]][["func"]]        <- c(rep(0, previousOutputLength), bn[[gene]][["func"]])
    bn[[gene]][["expr"]]        <- sumOfProducts(gene, bn[[gene]][["incoming"]], bn[[gene]][["func"]])
    return(bn)
}
#################
###### RND ######
#################
constRNDSelfLoop <- function(bn, gene){
    previousIncoming            <- bn[[gene]][["incoming"]]
    previousIncoming[1]         <- gene
    bn[[gene]][["incoming"]]    <- previousIncoming
    bn[[gene]][["expr"]]        <- sumOfProducts(gene, bn[[gene]][["incoming"]], bn[[gene]][["func"]])
    return(bn)
}
augmRNDSelfLoop <- function(bn, gene){
    previousIncoming            <- bn[[gene]][["incoming"]]
    previousIncoming            <- c(gene, previousIncoming)
    previousOutputLength        <- length(bn[[gene]][["func"]])
    bn[[gene]][["k"]]           <- bn[[gene]][["k"]] + 1
    bn[[gene]][["incoming"]]    <- previousIncoming
    bn[[gene]][["func"]]        <- c(bn[[gene]][["func"]], rndBooleanVector(previousOutputLength, bn[[gene]][["p"]]))
    bn[[gene]][["expr"]]        <- sumOfProducts(gene, bn[[gene]][["incoming"]], bn[[gene]][["func"]])
    return(bn)
}

addSelfLoops <- function(bn, numberOfSelfLoopsToAdd, SELF_LOOP_TYPE_FUNCTION){
    if (numberOfSelfLoops(bn) + numberOfSelfLoopsToAdd > length(bn)){
        stop("Impossible to add the required number of selfloops, its number is greater than the number of genes!")
    }
    i=1
    while(numberOfSelfLoopsToAdd > 0 & i <= length(bn)){
        if (!(bn[[i]][["id"]] %in% bn[[i]][["incoming"]])){ #non è un selfloop
            bn <- SELF_LOOP_TYPE_FUNCTION(bn, i)
            numberOfSelfLoopsToAdd <- numberOfSelfLoopsToAdd - 1
        }
        i = i + 1
    }
    return (bn)
}

# No. of selfloops in a Boolean network
getSelfLoops <- function(bn){
    sl <- sapply(bn, FUN= function(i) {if (i$id %in% i$incoming) return (i$id) else NA})
    sl <- sl[!is.na(sl)]
    return(sl)
}

# No. of selfloops in a Boolean network
numberOfSelfLoops <- function(bn){
    sl <- sapply(bn, FUN= function(i) i$id %in% i$incoming)
    return(sum(sl, na.rm = TRUE))
}

computePseudoAttractor <- function(att){
    means <- colMeans(att)
    return(sapply(means, FUN=function(x) ifelse(x>=0.5, 1, 0)))
}

#the set of nodes that always assume the same value—either 0 or 1—along all the attractors of the BN.
commonSea <- function(attractors){
    noAttractors <- length(attractors$attractors)

    df <- computePseudoAttractor(getAttractorSequence(attractors, 1))
    if (noAttractors > 1){
        for (i in seq(2,noAttractors)){
           df <- rbind(df, computePseudoAttractor(getAttractorSequence(attractors, i)))
        }  
        temp <- colSums(df, na.rm=FALSE)/nrow(df)
    } else {
        temp <- df
    }
    commonSeaOnes   <- which(temp==1)
    commonSeaZeros  <- which(temp==0)
    specificIsland <- temp[which(temp!=1 & temp!=0)]
    return(
            list("commonSeaSize"  = length(commonSeaOnes) + length(commonSeaZeros),
                 "commonSeaOnes"  = commonSeaOnes,
                 "commonSeaZeros" = commonSeaZeros,
                 "specificIslandSize" = length(specificIsland),
                 "specificIsland"= specificIsland
                )
            )
}

# Format conversion, from my description to BoolNet object, passing from filesystem
toBoolNet <- function(bn, filename) {
    saveNetworkToFile(bn, filename)
    return(loadNetwork(filename))
}
# Change temporarily the bit of boolean vector to respect the fixed nodes of the Boolean network
applyFixedGenes <- function(bn, booleanState) {
    for (i in seq_len(length(booleanState))){
        if ( bn$fixed[i] != -1){
            booleanState[i] <- bn$fixed[i]
        }
    }
    return(booleanState)
}

#Retrieve attractors starting from bagOfStates
simulateNet <- function(bn, initialStates) {
    initialStatesWithFixedGenesPerNetwork <- lapply(initialStates, FUN=function(x) applyFixedGenes(bn,x))
    atts <- getAttractors(bn, 
                            type = "synchronous", 
                            method = "chosen", 
                            startStates = initialStatesWithFixedGenesPerNetwork)
    return(atts)
}

#Functions that returns the Boolean functions with no irrelevant genes
booleanFunctionsWithNoIrrelevantGenes <- function(k){
    # Generate a list of all assignments of n variables with N possible values
    allcombn <- function(N,n)
    {
        rownum = N^n
        sapply(n:1,function(i)
            {
                    rep(seq_len(N),each=N^(i-1),len=rownum)
                })
    }
    table <- allcombn(2,k)-1
    numOfPossibleBoolFuns <- 2^(2^k)

    functionsWithNoIrrelevantGenes <- list()

    for (decimalNum in seq(0,numOfPossibleBoolFuns-1)){
        func <- fromIntegerToBitVector(integer=decimalNum, numOfBits=2^k)
        dropGenes <- apply(table,2,function(gene)
                  # determine all genes that have no influence on the results,
                  # i.e. the result column is equal for 0 values and 1 values
                          {
                            (identical(func[gene==1],
                                  func[gene==0]))
                          })
        if (sum(dropGenes) == 0){
            functionsWithNoIrrelevantGenes[[length(functionsWithNoIrrelevantGenes) + 1]] <- func
        } 
    }
    return (functionsWithNoIrrelevantGenes)
}


#bn <- rbn(5,3,0.9,selfLoops=FALSE)
#print("TOPOLOGY-PRE")
#print(lapply(bn, `[[`, "expr"))
#print("TOPOLOGY-POST")


#bn <- constANDSelfLoop(bn,5)

#bn <- augmANDSelfLoop(bn,1)
#print(numberOfSelfLoops(bn))

#bn <- addSelfLoops(bn, 3, constANDSelfLoop)
#print(numberOfSelfLoops(bn))

#print(lapply(bn, `[[`, "expr"))

#bb <- toBoolNet(bn, "prova/pippo")
#print(bb)
#plotNetworkWiring(bb)

#atts <- getAttractors(bb)
#print(atts$attractors)
#print(length(atts$attractors))
#commonSea(atts)


#sample(k2Functions,100,replace=TRUE)
#rbnSubset <- function(n, k, p, selfLoops=FALSE, setOfAllowedFunctions){
#sample(booleanFunctionsWithNoIrrelevantGenes(2),10,replace=TRUE)


#allk2Functions <- lapply(seq(0,15), FUN=fromIntegerToBitVector, numOfBits=4)



#booleanFunctionsWithNoIrrelevantGenes(2)
