
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
    print(booleanFunction)
    for (index in seq(1,length(booleanFunction))){
        if (booleanFunction[index] == 1){
            print(index)
            print(booleanFunction[index])
            minterm <- rev(as.numeric(intToBits(index-1))[1:numGenes])
            print(minterm)
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

rbn <- function(n, k){
    bn <- vector("list", length = n)
    nodes <- c(1:n)
    for (i in (seq(1,n))){
        bn[[i]][["k"]]          <- k
        bn[[i]][["incoming"]]   <- sample(nodes[!nodes == i], k, replace=FALSE)
        bn[[i]][["func"]]       <- rndBooleanState(2^bn[[i]][["k"]])
        bn[[i]][["expr"]]       <- sumOfProducts(i, bn[[i]][["incoming"]], bn[[i]][["func"]])
    }
    return(bn)
} 

rbn(10,3)
