library(BoolNet)
source("bOOLnET.R")

simulateNet <- function(bn, initialStates) {
    initialStatesWithFixedGenesPerNetwork <- lapply(initialStates, FUN=function(x) applyFixedGenes(bn,x))
    atts <- getAttractors(bn, 
                            type = "synchronous", 
                            method = "chosen", 
                            startStates = initialStatesWithFixedGenesPerNetwork)
    return(commonSea(atts))
}


set.seed(2)
k       <- 2
bias    <- 0.5
noNetworks <- 100

stats <- list()
for (noNodes in seq(5,14)){
    print(noNodes)
    res <- vector("integer", length = noNetworks)

    for(i in seq_len(noNetworks)){
        noInitialStates <- (2^noNodes)/2
        initialStates <- lapply(rep(1, noInitialStates), function(x) rndBooleanState(noNodes))

        net <- rbn(noNodes, k, bias, selfLoops=FALSE)
        #saveRDS(net, file = glue('{folder}/bn_{i}.RData'))
        bn  <- toBoolNet(net, glue('{subFolderSL_0}/bn_{i}'))

            
        res[i] <-   simulateNet(bn, initialStates)$commonSeaSize
       
    }
    stats[[length(stats) +1 ]] <- res
}

boxplot(stats)

#lapply(res, `[[`, "commonSeaSize")
#pippo<-readRDS("net_1.RData")
#print(addSelfLoops(pippo, 3, constORSelfLoop))

#saveNetwork(net, "prova/rete")