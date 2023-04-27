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

hier.clustering <- function(dist){
    res.hc <- hclust(d = dist, method = "complete")  
    return (res.hc)
    #return(fviz_dend(res.hc, cex = 0.5))
}

beforeAfterAnalysis <- function(df.complete, df.filtered, title, path, filename){
    dist.complete <- pegas::dist.hamming(df.complete)
    dist.filtered <- pegas::dist.hamming(df.filtered)
    
    # DENDROGRAMMA #
    pdf(glue('{path}dendro_{filename}'))
    dendro.plots <- list()  # new empty list
    dendro.plots[[1]] <- fviz_dend(hier.clustering(dist.complete), cex = 0.5, main=glue('complete - {title}'))
    dendro.plots[[2]] <- fviz_dend(hier.clustering(dist.filtered), cex = 0.5, main=glue('filtered - {title}'))
    multiplot(plotlist = dendro.plots, cols = 2)
    dev.off()  
    ################  

    # DISTRIBUZIONI DI DISTANZE #
    pdf(glue('{path}dist_{filename}'))
    distrib.plots <- list()  # new empty list
    distrib.complete <- unlist(dist.complete)
    distrib.plots[[1]] <- hist(distrib.complete, 
                    breaks=0:(max(distrib.complete)+1), 
                    xaxt="n", 
                    right=FALSE, 
                    freq=TRUE,
                    cex.lab=1, 
                    cex.axis=1.1,
                    cex.main=0.9,
                    ylab="Absolute Frequency",
                    xlab="Hamming distance",
                    main=glue('COMPLETE'))
        #tmp$counts=tmp$counts/sum(tmp$counts)
    axis(1, at=distrib.plots[[1]]$mids, labels=0:max(distrib.complete), cex.axis=1.1)

    distrib.filtered <- unlist(dist.filtered)
    distrib.plots[[2]] <- hist(distrib.filtered, 
                    breaks=0:(max(distrib.filtered)+1), 
                    xaxt="n", 
                    right=FALSE, 
                    freq=TRUE,
                    cex.lab=1, 
                    cex.axis=1.1,
                    cex.main=0.9,
                    ylab="Absolute Frequency",
                    xlab="Hamming distance",
                    main=glue('FILTERED'))
        #tmp$counts=tmp$counts/sum(tmp$counts)
    axis(1, at=distrib.plots[[2]]$mids, labels=0:max(distrib.filtered), cex.axis=1.1)
    multiplot(plotlist = distrib.plots, cols = 3)
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
for (SLNUMBER in c(1))#,2,3,4,5,10,20))
{
    for (SLTYPE in c("_augmAND","_augmOR")){
        numberOfAttractors <- c()
        for (NET in seq(1,1)){
            specific_path <- glue('{path}/sl{SLNUMBER}{SLTYPE}/atts/bn_{NET}_pseudoAtts.csv')
            #print(specific_path)
            df <- read.csv(specific_path, header=TRUE)

            condition <- seq_len(nrow(df))    
            for (COLUMN in seq(1,SLNUMBER)){
                if (SLTYPE == "_augmAND"){
                    condition <- intersect(condition, which(df[,COLUMN] == 0))
                } else {
                    condition <- intersect(condition, which(df[,COLUMN] == 1))
                }
            }
            #before <- nrow(df)
            no_filtered_attrs <- nrow(df[condition, ])
            #print(glue('BEFORE: {before}, AFTER {no_filtered_attrs}'))
            numberOfAttractors <- c(numberOfAttractors, no_filtered_attrs)

            if (no_filtered_attrs > 1 ){
                #saveDendrogram(df[condition, ], glue('dendro_{NET}_sl{SLNUMBER}{SLTYPE}_{type}.pdf'))
                beforeAfterAnalysis(df,df[condition, ],  glue('bn{NET}_sl{SLNUMBER}{SLTYPE}'), sub.path, glue('bn{NET}_sl{SLNUMBER}{SLTYPE}_{type}.pdf'))
            }

        }
        res[[length(res) + 1]] <- numberOfAttractors
        names(res)[length(res)] <- glue('{SLNUMBER}{SLTYPE}')
    }
}

#print(res)
#print(names(res))

pdf(glue('filtered_noPseudoAttrs_{type}.pdf'))
par( mar = c(6.5, 4, 2, 2))

boxplot(res, xaxt = "n",  ylab="no. of (pseudo)attractors")
axis(1, at = seq(1,length(res)), las = 2,labels = names(res)) # axis, ticks
dev.off()

#pdf(glue('DISTANCES_{type}_sl_0.pdf'))
#dev.off()
