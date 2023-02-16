library(BoolNet)
source("../bOOLnET.R")

set.seed(11)
k       <- 2
bias    <- 0.5
noNetworks <- 100

stats <- list()
stats_num_attractors <- list()
stats_atts_mean_size <- list()

start <- 50
end <- 200
mainFolder <- glue('testing_all_functions')
dir.create(mainFolder)
#delta <- 500
### CASO SENZA TRUE E FALSE ###
#allowedFunctions <- lapply(seq(0,15), FUN=fromIntegerToBitVector, numOfBits=4)
#allowedFunctions
#allowedFunctions[[1]] <- NULL
#allowedFunctions[[length(allowedFunctions)]] <- NULL
####
### CASO SENZA IRRELEVANT GENES ###
#allowedFunctions <- booleanFunctionsWithNoIrrelevantGenes(k=2)
###
nodesNumberVector <- c(50,100,200)
for (noNodes in nodesNumberVector){
    print(noNodes)
    res <- vector("integer", length = noNetworks)
    res_num_attractors <- vector("integer", length = noNetworks)
    res_mean_attrs_size <- vector("numeric", length = noNetworks)
    for(i in seq_len(noNetworks)){
        net <- rbn(n=noNodes, k=k, p=bias, selfLoops=FALSE)
        #net <- rbnSubset(n=noNodes, k=k, p=bias, selfLoops=FALSE, setOfAllowedFunctions=allowedFunctions)

        bn  <- toBoolNet(net, glue('{mainFolder}/bn_{noNodes}_{i}'))

        #bn <- generateRandomNKNetwork(n=noNodes,
         #                               k=k,
         #                               functionGeneration = "biased", 
		#								zeroBias = bias,
		#								readableFunctions=TRUE,
		#								noIrrelevantGenes=TRUE)

        atts <- getAttractors(bn, 
                            type = "synchronous", 
                            method="random",
                            startStates=1000) 
                            #startStates=as.integer((2^noNodes)/1000000)) 
        #atts <- getAttractors(bn, type = "synchronous") 
        res[i] <- (commonSea(atts)$commonSeaSize)/noNodes
        res_num_attractors[i] <- length(atts$attractors)
        res_mean_attrs_size[i] <- mean(sapply(atts$attractors, FUN=function(x) ncol(x$involvedStates)))
    }
    stats[[length(stats) +1 ]] <- res
    stats_num_attractors[[length(stats_num_attractors) +1]] <- res_num_attractors
    stats_atts_mean_size[[length(stats_atts_mean_size) +1]] <- res_mean_attrs_size
}

pdf(glue('NUOVI_RISULTATI_CommonSea_{start}-{end}_100_all_functions.pdf'))
boxplot(stats,xaxt = "n", ylab="Common sea size / no. nodes", main="Size progression of the common sea (all functions allowed)")
#axis(1, at = seq(1,((end-start)/delta)+1), labels = seq(start, end,by=delta)) # axis, ticks
axis(1, at = seq(1,length(nodesNumberVector)), labels = nodesNumberVector) # axis, ticks
mtext('no. nodes', side=1, line=3)
dev.off()

pdf(glue('NUOVI_RISULTATI_CommonSea_{start}-{end}_100_all_functions_attractors.pdf'))
boxplot(stats_num_attractors,log="y",xaxt = "n", ylab="No. of attractors", main="No. of attractors")
#axis(1, at = seq(1,((end-start)/delta)+1), labels = seq(start, end,by=delta)) # axis, ticks
axis(1, at = seq(1,length(nodesNumberVector)), labels = nodesNumberVector) # axis, ticks
mtext('no. nodes', side=1, line=3)
dev.off()


pdf(glue('NUOVI_RISULTATI_CommonSea_{start}-{end}_100_all_functions_attractors_period.pdf'))
boxplot(stats_atts_mean_size, log="y",xaxt = "n", ylab="Mean of the period", main="Mean of the attractors' period")
#axis(1, at = seq(1,((end-start)/delta)+1), labels = seq(start, end,by=delta)) # axis, ticks
axis(1, at = seq(1,length(nodesNumberVector)), labels = nodesNumberVector) # axis, ticks
mtext('no. nodes', side=1, line=3)
dev.off()

#lapply(res, `[[`, "commonSeaSize")
#pippo<-readRDS("net_1.RData")
#print(addSelfLoops(pippo, 3, constORSelfLoop))

#saveNetwork(net, "prova/rete")