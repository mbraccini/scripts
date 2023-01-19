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

get_fun <- function(fun){
             fun <- deparse(fun)
             chunk <- tail(fun, 1)
             words <- strsplit(chunk, "\"")[[1]]
             return(words[2])
           }

set.seed(2)
k       <- 2
bias    <- 0.5
noInitialStates <- 10
noNodes <- 10
mainFolder <- glue('n{noNodes}k{k}p{gsub(".","",bias,fixed=TRUE)}')
dir.create(mainFolder)
initialStates <- lapply(rep(1, noInitialStates), function(x) rndBooleanState(noNodes))
noNetworks <- 1
subFolderSL_0 <- glue('{mainFolder}/sl0')
dir.create(subFolderSL_0)

res <- vector("list", length = noNetworks)
for(i in seq_len(noNetworks)){
    net <- rbn(noNodes, k, bias, selfLoops=FALSE)
    #saveRDS(net, file = glue('{folder}/bn_{i}.RData'))
    bn  <- toBoolNet(net, glue('{subFolderSL_0}/bn_{i}'))

    #### SIMULATE NET ####
    simulateNet(bn, initialStates)
    ######################

    for (noSelfLoops in c(2, 4, 6)){
        toAdd <- noSelfLoops - numberOfSelfLoops(net)

        for (selfLoopType in c(constANDSelfLoop)){
            slType          <- as.character(substitute(selfLoopType))
            print(get_fun(slType))
            subFolderSL_n    <- glue('{mainFolder}/sl{noSelfLoops}_{deparse(substitute(selfLoopType))}')
            dir.create(subFolderSL_n)

            net <- addSelfLoops(net, toAdd, selfLoopType)
            bn  <- toBoolNet(net, glue('{subFolderSL_n}/bn_{i}'))

            #### SIMULATE NET ####
            simulateNet(bn, initialStates)
            ######################
         
        }
    }
}
res

#lapply(res, `[[`, "commonSeaSize")
#pippo<-readRDS("net_1.RData")
#print(addSelfLoops(pippo, 3, constORSelfLoop))

#saveNetwork(net, "prova/rete")