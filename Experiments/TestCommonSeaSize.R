library(BoolNet)
source("../bOOLnET.R")



set.seed(42)
k       <- 2
bias    <- 0.5
noNetworks <- 100

stats <- list()
start <- 30
end <- 40
mainFolder <- glue('testing')
dir.create(mainFolder)

for (noNodes in seq(start,end,by=5)){
    print(noNodes)
    res <- vector("integer", length = noNetworks)

    for(i in seq_len(noNetworks)){
        net <- rbn(n=noNodes, k=k, p=bias, selfLoops=FALSE)
        bn  <- toBoolNet(net, glue('{mainFolder}/bn_{noNodes}_{i}'))

        #bn <- generateRandomNKNetwork(n=noNodes, k=k)
        atts <- getAttractors(bn, type = "synchronous", method="random",startStates=as.integer((2^noNodes)/1000000)) 
        #atts <- getAttractors(bn, type = "synchronous") 
        res[i] <- (commonSea(atts)$commonSeaSize)/noNodes
    }
    stats[[length(stats) +1 ]] <- res
}

pdf(glue('CommonSeaSizeProgression{start}-{end}_1000000_senzaatuoanelli.pdf'))
boxplot(stats,xaxt = "n", ylab="Common sea size / no. nodes", main="Size progression of the common sea")
axis(1, at = seq(1,((end-start)/5)+1), labels = seq(start, end,by=5)) # axis, ticks
mtext('no. nodes', side=1, line=3)
dev.off()

#lapply(res, `[[`, "commonSeaSize")
#pippo<-readRDS("net_1.RData")
#print(addSelfLoops(pippo, 3, constORSelfLoop))

#saveNetwork(net, "prova/rete")