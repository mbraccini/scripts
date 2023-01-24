library(glue)

path="n100k2p05/"

res     <- list()
labels  <- list()
sl0     <- read.table(glue('{path}commonSeaSize_sl_0.txt'))[,1]
res[[1]]<- sl0
labels[[1]]  <- "0"
for (SLNUMBER in c(5,10,20))
{
    res[[length(res)+1]] <- read.table(glue('{path}commonSeaSize_sl_{SLNUMBER}_augmAND.txt'))[,1]
    labels[[length(labels)+1]] <- glue('{SLNUMBER}_And')
    res[[length(res)+1]] <- read.table(glue('{path}commonSeaSize_sl_{SLNUMBER}_augmOR.txt'))[,1]
    labels[[length(labels)+1]] <- glue('{SLNUMBER}_Or')
    res[[length(res)+1]] <- read.table(glue('{path}commonSeaSize_sl_{SLNUMBER}_augmRND.txt'))[,1]
    labels[[length(labels)+1]] <- glue('{SLNUMBER}_Rnd')
}

boxplot(res, xaxt = "n", main="n=100, k=2, p=0.5, samples=100, initialStates=1000,\nselfLoopsType= augm{AND,OR,RND}", ylab="Common Sea Size (no. of nodes)")
axis(1, at = seq(1,length(res)), las = 2,labels = labels) # axis, ticks
dev.off()