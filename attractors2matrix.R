#Function - from attractors and gene numbers retrieve a list of matrices with columns that are the states that compose the attractors of the net
getMatricesAttractors <- function(attractors, numGens){
    matrixAttractorList <- list()
    #print(attractors$attractors[[1]]$involvedStates) #A matrix containing the states that make up the attractor. Each column represents one state.  The entries are decimal numbers that internally represent the states of the genes
    for (i in c(1:length(attractors$attractors))){ #per il numero di attrattori
        matrixAttractorList[[i]] <- attractors$attractors[[i]]$involvedStates
    }
    attsNames <- attractors$stateInfo$genes  #ci appiccichiamo anche il nome dei geni
    intToBitVector <- function(attractorColumn, numBit){ 
        v <- as.numeric(intToBits(attractorColumn))[1:numBit]
        names(v) <- attsNames
        return (v)
    }
    funForAMatrixDescribingTheAttractor <- function(attractorsMatrix,numBits){ apply(attractorsMatrix,2,intToBitVector, numBit=numBits)} #2 significa applicalo sulle colonne
    syncAttractors <- lapply(matrixAttractorList,funForAMatrixDescribingTheAttractor,numBits=numGens) 
    return (syncAttractors)
}