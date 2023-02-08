library(BoolNet)
source("bOOLnET.R")



set.seed(42)
k       <- 2
bias    <- 0.5
noNetworks <- 100

stats <- list()
start <- 20
end <- 100
mainFolder <- glue('testing')
dir.create(mainFolder)
delta <- 5
for (noNodes in seq(start,end,by=delta)){
    print(noNodes)
    res <- vector("integer", length = noNetworks)

    for(i in seq_len(noNetworks)){
        #net <- rbn(n=noNodes, k=k, p=bias, selfLoops=FALSE)
        #bn  <- toBoolNet(net, glue('{mainFolder}/bn_{noNodes}_{i}'))

        bn <- generateRandomNKNetwork(n=noNodes,
                                        k=k,
                                        functionGeneration = "biased", 
										zeroBias = bias,
										readableFunctions=TRUE,
										noIrrelevantGenes=TRUE)
        atts <- getAttractors(bn, 
                            type = "synchronous", 
                            method="random",
                            startStates=1000) 
                            #startStates=as.integer((2^noNodes)/1000000)) 
        #atts <- getAttractors(bn, type = "synchronous") 
        res[i] <- (commonSea(atts)$commonSeaSize)/noNodes
    }
    stats[[length(stats) +1 ]] <- res
}

pdf(glue('CommonSeaSizeProgression{start}-{end}_1000_noIrrelevantGenes_TRUE_Boolnet_FRACTION.pdf'))
boxplot(stats,xaxt = "n", ylab="Common sea size / no. nodes", main="Size progression of the common sea")
axis(1, at = seq(1,((end-start)/delta)+1), labels = seq(start, end,by=delta)) # axis, ticks
mtext('no. nodes', side=1, line=3)
dev.off()

#lapply(res, `[[`, "commonSeaSize")
#pippo<-readRDS("net_1.RData")
#print(addSelfLoops(pippo, 3, constORSelfLoop))

#saveNetwork(net, "prova/rete")