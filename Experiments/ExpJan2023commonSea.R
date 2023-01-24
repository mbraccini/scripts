library(BoolNet)
source("bOOLnET.R")


#Check if selfloops are in the commonSea passed (COmposed of all 1's or all 0's)
retrieveSelfLoopsInCommonSea <- function(myNet, commonSea){
        tempSelfLoopsNames  <- paste0("Gene",getSelfLoops(myNet))
        commonSeaONES       <- names(commonSea$commonSeaOnes)
        CommonSeaZEROS      <- names(commonSea$commonSeaZeros)
        inONES      <- tempSelfLoopsNames[tempSelfLoopsNames %in% commonSeaONES]
        inZEROS     <- tempSelfLoopsNames[tempSelfLoopsNames %in% CommonSeaZEROS]
        return(list("inONES"=inONES,"inZEROS"=inZEROS))
}

set.seed(5)
k       <- 2
bias    <- 0.5
noInitialStates <- 1000
noNodes <- 100
mainFolder <- glue('n{noNodes}k{k}p{gsub(".","",bias,fixed=TRUE)}')
dir.create(mainFolder)
initialStates <- lapply(rep(1, noInitialStates), function(x) rndBooleanState(noNodes))
noNetworks <- 100
subFolderSL_0 <- glue('{mainFolder}/sl0')
dir.create(subFolderSL_0)

for(i in seq_len(noNetworks)){
    net <- rbn(noNodes, k, bias, selfLoops=FALSE)
    #saveRDS(net, file = glue('{folder}/bn_{i}.RData'))
    bn  <- toBoolNet(net, glue('{subFolderSL_0}/bn_{i}'))

    #### SIMULATE NET ####
    atts <- simulateNet(bn, initialStates)
    #saveAttractors(atts, glue('bn_{i}'),glue('{subFolderSL_0}') )
    CS <- commonSea(atts)
    ######################
    write(CS$commonSeaSize,  file=glue('{mainFolder}/commonSeaSize_sl_0.txt'),  append=TRUE)

    for (noSelfLoops in c(5, 10, 20)){
        for (selfLoopType in c("OR","AND","RND")){
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
            res <- retrieveSelfLoopsInCommonSea(temp_net, CS)
            res_elab <- lapply(res, FUN=length)
            write(res_elab$inONES,  file=glue('{mainFolder}/slInOnes_{noSelfLoops}_{SL_TYPE_STRING}.txt'),  append=TRUE)
            write(res_elab$inZEROS, file=glue('{mainFolder}/slInZeros_{noSelfLoops}_{SL_TYPE_STRING}.txt'), append=TRUE)

            write(CS$commonSeaSize,  file=glue('{mainFolder}/commonSeaSize_sl_{noSelfLoops}_{SL_TYPE_STRING}.txt'),  append=TRUE)

        }
    }
}


#lapply(res, `[[`, "commonSeaSize")
#pippo<-readRDS("net_1.RData")
#print(addSelfLoops(pippo, 3, constORSelfLoop))

#saveNetwork(net, "prova/rete")