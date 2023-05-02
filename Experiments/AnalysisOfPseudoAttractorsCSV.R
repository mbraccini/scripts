library(glue)
library(pegas)
library("factoextra")
library(ggplot2)

#type = "allFunctions" 
type = "noTRUE_FALSE_XOR_XNOR"
path=glue('n100k2p05_pseudoAttractors_{type}_26Aprile23')
sub.path = glue('{path}/analysis/')
dir.create(sub.path)

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

res <- list()
#names(a)[1] <- c("biill")
for (SLNUMBER in c(1,2,3,4,5,10,20))
{
    for (SLTYPE in c("_augmAND","_augmOR")){
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
            before <- nrow(df)
            no_filtered_attrs <- nrow(df[condition, ])
            #print(glue('BEFORE: {before}, AFTER {no_filtered_attrs}'))
            numberOfAttractorsBEFORE    <- c(numberOfAttractorsBEFORE, before)
            numberOfAttractors          <- c(numberOfAttractors, no_filtered_attrs)

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
        
    }
}

#print(res)
#print(names(res))

pdf(glue('filtered_noPseudoAttrs_{type}.pdf'))
par( mar = c(8.5, 4, 2, 2))

boxplot(res , log="y", xaxt = "n",  ylab="no. of (pseudo)attractors + 1")
axis(1, at = seq(1,length(res)), las = 2,labels = names(res)) # axis, ticks
dev.off()

#pdf(glue('DISTANCES_{type}_sl_0.pdf'))
#dev.off()
