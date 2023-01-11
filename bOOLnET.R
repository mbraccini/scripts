
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


rbn <- function(n, k){
    bn <- vector("list", length = n)
    nodes <- c(1:n)
    for (i in (seq(1,n))){
        bn[[i]][["k"]] <- k
        bn[[i]][["incoming"]] <- sample(nodes[!nodes == i], k, replace=FALSE)
        bn[[i]][["func"]] <- rndBooleanState(2^bn[[i]][["k"]])
    }
    return(bn)
} 

rbn(10,3)