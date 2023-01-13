
library(glue)
#number of nodes in a network loaded with the load() function
numberOfGenes <- function(net){
    return (length(net[[1]])) # We read from the imported network its number of nodes
}

saveAttractors <- function(attractors, path="./"){
    for (i in seq(1,length(attractors$attractors))){
        att <- getAttractorSequence(attractors, i)
        write.table(att, glue('{path}att_{i}.csv'), sep = ",", col.names = colnames(att),row.names=FALSE)
    }
}


#Generate a RANDOM Boolean state of a specified length
rndBooleanState <- function(vectorLength){ return (sample(c(0,1), vectorLength, replace=TRUE)) }

sumOfProducts <- function(gene, incomingGenes, booleanFunction){
    numGenes <- length(incomingGenes)
    string <- glue('Gene{gene}, ')
    #print(booleanFunction)
    for (index in seq(1,length(booleanFunction))){
        if (booleanFunction[index] == 1){
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
    string <- gsub('.{3}$', '', string)
    return(string)
}

rbn <- function(n, k, selfLoops=FALSE){
    bn <- vector("list", length = n)
    nodes <- c(1:n)
    for (i in (seq(1,n))){
        bn[[i]][["k"]]          <- k
        if (selfLoops) {
            bn[[i]][["incoming"]]   <- sample(nodes, k, replace=FALSE)
        } else {
            bn[[i]][["incoming"]]   <- sample(nodes[!nodes == i], k, replace=FALSE)
        }
        bn[[i]][["func"]]       <- rndBooleanState(2^bn[[i]][["k"]])
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
constOrSelfLoop <- function(bn, gene){
    previousIncoming            <- bn[[gene]][["incoming"]]
    previousIncoming[1]         <- gene
    bn[[gene]][["incoming"]]    <- previousIncoming
    bn[[gene]][["func"]]        <- c(0 ,rep(1, 2^bn[[gene]][["k"]] - 1))
    bn[[gene]][["expr"]]        <- sumOfProducts(gene, bn[[gene]][["incoming"]], bn[[gene]][["func"]])
    return(bn)
}
augmOrSelfLoop <- function(bn, gene){
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
constAndSelfLoop <- function(bn, gene){
    previousIncoming            <- bn[[gene]][["incoming"]]
    previousIncoming[1]         <- gene
    bn[[gene]][["incoming"]]    <- previousIncoming
    bn[[gene]][["func"]]        <- c(rep(0, 2^bn[[gene]][["k"]] - 1), 1)
    bn[[gene]][["expr"]]        <- sumOfProducts(gene, bn[[gene]][["incoming"]], bn[[gene]][["func"]])
    return(bn)
}
augmAndSelfLoop <- function(bn, gene){
    previousIncoming            <- bn[[gene]][["incoming"]]
    previousIncoming            <- c(gene, previousIncoming)
    previousOutputLength        <- length(bn[[gene]][["func"]])
    bn[[gene]][["k"]]           <- bn[[gene]][["k"]] + 1
    bn[[gene]][["incoming"]]    <- previousIncoming
    bn[[gene]][["func"]]        <- c(rep(0, previousOutputLength), bn[[gene]][["func"]])
    bn[[gene]][["expr"]]        <- sumOfProducts(gene, bn[[gene]][["incoming"]], bn[[gene]][["func"]])
    return(bn)
}



bn <- rbn(10,4, selfLoops=FALSE)
print("TOPOLOGY-PRE")
print(lapply(bn, `[[`, "expr"))
print("TOPOLOGY-POST")
print(bn[[10]])

bn <- augmAndSelfLoop(bn,10)
print(bn[[10]])
print(lapply(bn, `[[`, "expr"))
saveNetworkToFile(bn, "prova/pippo")
bb <- loadNetwork("prova/pippo")
print(bb)
plotNetworkWiring(bb)