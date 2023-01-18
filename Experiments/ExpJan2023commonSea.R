library(BoolNet)
source("bOOLnET.R")

set.seed(2)
k       <- 2
bias    <- 0.5
noInitialStates <- 1000
nodesList <- c(10,30,50)

results <- list()
for(noNodes in nodesList){
    folder <- glue('nodes_{noNodes}_sl_0')
    dir.create(folder)
    initialStates <- lapply(rep(1, noInitialStates), function(x) rndBooleanState(noNodes))
    ##########
    noNetworks <- 100
    res <- vector("integer", length = noNetworks)
    for(i in seq_len(noNetworks)){
        net <- rbn(noNodes, k, bias, selfLoops=FALSE)
        saveRDS(net, file = glue('{folder}/bn_{i}.RData'))
        bn  <- toBoolNet(net, glue('{folder}/bn_{i}'))

        initialStatesWithFixedGenesPerNetwork <- lapply(initialStates, FUN=function(x) applyFixedGenes(bn,x))
        atts <- getAttractors(bn, 
                            type = "synchronous", 
                            method = "chosen", 
                            startStates = initialStatesWithFixedGenesPerNetwork)
        res[i] <- commonSea(atts)$commonSeaSize
        
    }
    results[[length(results) + 1]] <- res
}

boxplot(results)
#lapply(res, `[[`, "commonSeaSize")
#pippo<-readRDS("net_1.RData")
#print(addSelfLoops(pippo, 3, constORSelfLoop))

#saveNetwork(net, "prova/rete")