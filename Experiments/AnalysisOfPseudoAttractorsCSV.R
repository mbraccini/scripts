library(glue)
library(pegas)
library("factoextra")
library(ggplot2)

#type = "allFunctions" 
type = "noTRUE_FALSE_XOR_XNOR"
path=glue('n100k2p05_pseudoAttractors_{type}_26Aprile23')
sub.path = glue('{path}/analysis/')
dir.create(sub.path)

commonSeaFromDataFrame <- function(df){
    
    temp <- colSums(df, na.rm=FALSE)/nrow(df)

    commonSeaOnes   <- which(temp==1)
    commonSeaZeros  <- which(temp==0)
    specificIsland <- temp[which(temp!=1 & temp!=0)]
    return(
            list("commonSeaSize"  = length(commonSeaOnes) + length(commonSeaZeros),
                 "commonSeaOnes"  = commonSeaOnes,
                 "commonSeaZeros" = commonSeaZeros,
                 "specificIslandSize" = length(specificIsland),
                 "specificIsland"= specificIsland
                )
            )
}

################################################
################################################
################################################
################################################
################################################
################################################
options(width=60)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
################################################
################################################
################################################
################################################
################################################
################################################


dendrogramma <- function(df.complete, df.filtered, title, path, filename){
    dist.complete <- pegas::dist.hamming(df.complete)
    dist.filtered <- pegas::dist.hamming(df.filtered)
    
    hc.complete <- hclust(d = dist.complete, method = "complete")  
    hc.filtered <- hclust(d = dist.filtered, method = "complete")  

    # DENDROGRAMMA #
    pdf(glue('{path}dendro_{filename}'))
    dendro.plots <- list()  # new empty list
    dendro.plots[[1]] <- fviz_dend(hc.complete, cex = 0.5, main=glue('COMPLETE - {title}'))
    dendro.plots[[2]] <- fviz_dend(hc.filtered, cex = 0.5, main=glue('FILTERED - {title}'))
    multiplot(plotlist = dendro.plots, cols = 2)
    dev.off()  
    ################  
    
}

generateDistributionsOfDistancesPlots <- function(distrib.complete, distrib.filtered, title, path, filename){
    # DISTRIBUZIONI DI DISTANZE #
        pdf(glue('{path}dist_{filename}'))
        par(mfrow=c(1,2))
        t.complete <- hist(distrib.complete, 
                        breaks=0:(max(distrib.complete)+1), 
                        xaxt="n", 
                        right=FALSE, 
                        freq=TRUE,
                        cex.lab=1, 
                        cex.axis=1.1,
                        cex.main=0.9,
                        ylab="Absolute Frequency",
                        xlab="Hamming distance",
                        main=glue('COMPLETE: {title}'))
            #tmp$counts=tmp$counts/sum(tmp$counts)
        axis(1, at=t.complete$mids, labels=0:max(distrib.complete), cex.axis=1.1)
        t.filtered <- hist(distrib.filtered, 
                        breaks=0:(max(distrib.filtered)+1), 
                        xaxt="n", 
                        right=FALSE, 
                        freq=TRUE,
                        cex.lab=1, 
                        cex.axis=1.1,
                        cex.main=0.9,
                        ylab="Absolute Frequency",
                        xlab="Hamming distance",
                        main=glue('FILTERED: {title}'))
            #tmp$counts=tmp$counts/sum(tmp$counts)
        axis(1, at=t.filtered$mids, labels=0:max(distrib.filtered), cex.axis=1.1)
        dev.off()  
}

#saveDendrogram <- function(df, filename){
#        res.dist <- pegas::dist.hamming(df)
#        print(hist(unlist(res.dist)))
#        res.hc <- hclust(d = res.dist, method = "complete")  
#        pdf(filename)
#        p <-fviz_dend(res.hc, cex = 0.5)
        #print(p)
        #print(p)
#        plots <- list()  # new empty list
#        plots[[1]] <- p
#        plots[[2]] <- p
#        multiplot(plotlist = plots, cols = 2)

#        dev.off()
#}

res             <- list()
resCommonSea    <- list()
#names(a)[1] <- c("biill")

#0 autoanelli
numberOfAttractorsBEFORE    <- c()
commonSeaBEFORE             <- c()
for (NET in seq(1,100)){
    path_0sl <- glue('{path}/sl0/atts/bn_{NET}_pseudoAtts.csv')
    df <- read.csv(path_0sl, header=TRUE)

    beforeNumAtts <- nrow(df)
    numberOfAttractorsBEFORE    <- c(numberOfAttractorsBEFORE, beforeNumAtts)

    if(beforeNumAtts > 1){
        beforeCommonSea <- commonSeaFromDataFrame(df)$commonSeaSize
        commonSeaBEFORE    <- c(commonSeaBEFORE, beforeCommonSea)
    }
}
res[[length(res) + 1]] <- numberOfAttractorsBEFORE + 1
names(res)[length(res)] <- glue('sl0')

resCommonSea[[length(resCommonSea) + 1]] <- commonSeaBEFORE 
names(resCommonSea)[length(resCommonSea)] <- glue('sl0')

for (SLNUMBER in c(1,2,3,4,5,10,20))
{
    for (SLTYPE in c("_augmAND","_augmOR")){
        commonSeaBEFORE             <- c()
        commonSea                   <- c()
        numberOfAttractorsBEFORE    <- c()
        numberOfAttractors          <- c()
        distrib.complete            <- c()
        distrib.filtered            <- c()
        for (NET in seq(1,100)){
            
            specific_path <- glue('{path}/sl{SLNUMBER}{SLTYPE}/atts/bn_{NET}_pseudoAtts.csv')
            df <- read.csv(specific_path, header=TRUE)

            net.path = glue('{path}/analysis/bn_{NET}/')
            dir.create(net.path)

            condition <- seq_len(nrow(df))    
            for (COLUMN in seq(1,SLNUMBER)){
                if (SLTYPE == "_augmAND"){
                    condition <- intersect(condition, which(df[,COLUMN] == 0))
                } else {
                    condition <- intersect(condition, which(df[,COLUMN] == 1))
                }
            }
            #NUMBER OF PSEUDOTTRACTORS
            beforeAtts  <- nrow(df)
            after   <- df[condition, ]
            no_filtered_attrs <- nrow(after)
            #print(glue('BEFORE: {before}, AFTER {no_filtered_attrs}'))
            numberOfAttractorsBEFORE    <- c(numberOfAttractorsBEFORE, beforeAtts)
            numberOfAttractors          <- c(numberOfAttractors, no_filtered_attrs)
            
            #COMMON SEA
            if(beforeAtts > 1 && no_filtered_attrs > 1){
                #print("un solo pseudoattrattore")
                beforeCommonSea     <- commonSeaFromDataFrame(df)$commonSeaSize
                commonSeaBEFORE     <- c(commonSeaBEFORE, beforeCommonSea)
            #}
            #if(no_filtered_attrs > 1){
                commonSeaAFTER       <- commonSeaFromDataFrame(after)$commonSeaSize
                commonSea           <- c(commonSea, commonSeaAFTER)
            }
            #if (no_filtered_attrs > 1 ){
            #    dendrogramma(df,
            #                df[condition, ],  
            #                glue('bn{NET}_sl{SLNUMBER}{SLTYPE}'), 
            #                net.path, 
            #                glue('bn{NET}_sl{SLNUMBER}{SLTYPE}_{type}.pdf'))
            #    distrib.complete <- c(distrib.complete, pegas::dist.hamming(df))
            #    distrib.filtered <- c(distrib.filtered, pegas::dist.hamming(df[condition, ]))
            #}

        }
        #generateDistributionsOfDistancesPlots(distrib.complete, 
        #                                    distrib.filtered, 
        #                                    glue('sl{SLNUMBER}{SLTYPE}'), 
        #                                    sub.path, 
        #                                    glue('sl{SLNUMBER}{SLTYPE}_{type}.pdf'))

        
        res[[length(res) + 1]] <- numberOfAttractorsBEFORE + 1
        names(res)[length(res)] <- glue('before_{SLNUMBER}{SLTYPE}')
        res[[length(res) + 1]] <- numberOfAttractors + 1
        names(res)[length(res)] <- glue('{SLNUMBER}{SLTYPE}')
        

        resCommonSea[[length(resCommonSea) + 1]] <- commonSeaBEFORE 
        names(resCommonSea)[length(resCommonSea)] <- glue('before_{SLNUMBER}{SLTYPE}')
        resCommonSea[[length(resCommonSea) + 1]] <- commonSea 
        names(resCommonSea)[length(resCommonSea)] <- glue('{SLNUMBER}{SLTYPE}')
    }
}

#print(res)
#print(names(res))

pdf(glue('filtered_noPseudoAttrs_{type}.pdf'))
par( mar = c(9, 4, 2, 2))

boxplot(res , log="y", xaxt = "n",  ylab="no. of (pseudo)attractors + 1")
axis(1, at = seq(1,length(res)), las = 2,labels = names(res)) # axis, ticks
dev.off()



pdf(glue('filtered_commonSea_{type}_greater1_both.pdf'))
par( mar = c(9, 4, 4, 2))

boxplot(resCommonSea , xaxt = "n",  ylab="Common Sea size",ylim=c(0,100),main="Considering RBNs with no. of pseudoattractors greater \n than 1 in both conditions (before and after filtering)")
axis(1, at = seq(1,length(resCommonSea)), las = 2,labels = names(resCommonSea)) # axis, ticks
dev.off()
#pdf(glue('DISTANCES_{type}_sl_0.pdf'))
#dev.off()
