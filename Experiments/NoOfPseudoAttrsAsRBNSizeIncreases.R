#!/usr/bin/env Rscript
library(BoolNet)
source("../bOOLnET.R")
args = commandArgs(trailingOnly=TRUE)

# Esempio: nohup Rscript ExpJan2023commonSea.R TRUE 876 &
ALL_FUNCTIONS <- as.logical(args[1])
seed = as.integer(args[2])
set.seed(seed)
print(glue('seed: {seed}'))
k       <- 2
bias    <- 0.5
noInitialStates <- 1000
noNetworks <- 100


if (!ALL_FUNCTIONS){
    ### RBN Subset ###
    allowedFunctions <- lapply(seq(0,15), FUN=fromIntegerToBitVector, numOfBits=4)
    allowedFunctions[[7]] <- NULL  #XOR
    allowedFunctions[[9]] <- NULL   #XNOR
    allowedFunctions[[1]] <- NULL   #FALSE
    allowedFunctions[[length(allowedFunctions)]] <- NULL    #TRUE
    print("allowed functions")
    allowedFunctions
    #################
    expString <- "_noTRUE_FALSE_XOR_XNOR"
} else {
    print("all functions")
    expString <- "_allFunctions"
}


for (noNodes in c(300, 500)){
    mainFolder <- glue('n{noNodes}k{k}p{gsub(".","",bias,fixed=TRUE)}_{expString}_noPseudoAttsAsSizeIncreases_3Aprile')
    dir.create(mainFolder)
    subFolderSL_0 <- glue('{mainFolder}/sl0')
    dir.create(subFolderSL_0)
    initialStates <- lapply(rep(1, noInitialStates), function(x) rndBooleanState(noNodes))


    for(i in seq_len(noNetworks)){
        if (!ALL_FUNCTIONS){
            net <- rbnSubset(n=noNodes, k=k, p=bias, selfLoops=FALSE, setOfAllowedFunctions=allowedFunctions)
        } else {
            net <- rbn(noNodes, k, bias, selfLoops=FALSE)
        }

        #saveRDS(net, file = glue('{folder}/bn_{i}.RData'))
        bn  <- toBoolNet(net, glue('{subFolderSL_0}/bn_{i}'))

        for (noSelfLoops in c(20)){
            for (selfLoopType in c("AND","RND")){
                result = switch(  
                selfLoopType,  
                "OR"= {
                    SL_FUN = augmORSelfLoop
                    SL_TYPE_STRING = "augmOR"
                },  
                "AND"= {
                    SL_FUN = augmANDSelfLoop
                    SL_TYPE_STRING = "augmAND"
                },
                "RND"= {
                    SL_FUN = augmRNDSelfLoop
                    SL_TYPE_STRING = "augmRND"
                },   
                stop(paste0("No handler for ", selfLoopType))
                )  

                subFolderSL_n    <- glue('{mainFolder}/sl{noSelfLoops}_{SL_TYPE_STRING}')
                dir.create(subFolderSL_n)

                temp_net <- addSelfLoops(net, noSelfLoops, SL_FUN)
                temp_bn  <- toBoolNet(temp_net, glue('{subFolderSL_n}/bn_{i}'))

                #### SIMULATE NET ####
                atts <- simulateNet(temp_bn, initialStates)
                #saveAttractors(atts, glue('bn_{i}'),glue('{subFolderSL_n}') )
                CS <- commonSea(atts)
                ######################
                distances <- distanceDistributionBetweenPseudoAttractors(atts)
                if (!is.null(distances)){
                    write.table(as.matrix(t(distances)), glue('{mainFolder}/distances_sl_{noSelfLoops}_{SL_TYPE_STRING}.csv'),sep=",",append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)
                    #per leggerlo basta usare read.table(file,sep=",",fill=TRUE)
                } 
                res <- retrieveSelfLoopsInCommonSea(temp_net, CS)
                res_elab <- lapply(res, FUN=length)

                noAttractors        <- length(atts$attractors)
                noPseudoAttractors  <- numberOfPseudoAttractors(atts)
                write(noAttractors,  file=glue('{mainFolder}/no_attractors_{noSelfLoops}_{SL_TYPE_STRING}.txt'),  append=TRUE)
                write(noPseudoAttractors,  file=glue('{mainFolder}/no_PseudoAttractors_{noSelfLoops}_{SL_TYPE_STRING}.txt'),  append=TRUE)

                write(res_elab$inONES,  file=glue('{mainFolder}/slInOnes_{noSelfLoops}_{SL_TYPE_STRING}.txt'),  append=TRUE)
                write(res_elab$inZEROS, file=glue('{mainFolder}/slInZeros_{noSelfLoops}_{SL_TYPE_STRING}.txt'), append=TRUE)

                write(CS$commonSeaSize,  file=glue('{mainFolder}/commonSeaSize_sl_{noSelfLoops}_{SL_TYPE_STRING}.txt'),  append=TRUE)
                write(length(CS$commonSeaOnes),  file=glue('{mainFolder}/commonSeaOnes_sl_{noSelfLoops}_{SL_TYPE_STRING}.txt'),  append=TRUE)

            }
        }
    }
}




