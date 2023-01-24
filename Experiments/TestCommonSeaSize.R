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
start <- 30
end <- 42
for (noNodes in seq(start,end,by=2)){
    print(noNodes)
    res <- vector("integer", length = noNetworks)

    for(i in seq_len(noNetworks)){
        bn <- generateRandomNKNetwork(n=noNodes, k=k, functionGeneration = "biased", zeroBias=bias)
        #atts <- getAttractors(bn, type = "synchronous") 
        atts <- getAttractors(bn, type = "synchronous", method="random",startStates=as.integer((2^noNodes) / 1000000) ) 
        res[i] <- (commonSea(atts)$commonSeaSize)/noNodes
    }
    stats[[length(stats) +1 ]] <- res
}

pdf(glue('CommonSeaSizeProgression{start}-{end}_sampled.pdf'))
boxplot(stats,xaxt = "n", ylab="Common sea size / no. nodes", main="Size progression of the common sea (1/1000000)")
axis(1, at = seq(1,((end-start)/2)+1), labels = seq(start, end,by=2)) # axis, ticks
mtext('no. nodes', side=1, line=3)
dev.off()

#lapply(res, `[[`, "commonSeaSize")
#pippo<-readRDS("net_1.RData")
#print(addSelfLoops(pippo, 3, constORSelfLoop))

#saveNetwork(net, "prova/rete")