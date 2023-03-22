library(glue)

path="n100k2p05_pseudoAttractors_noTRUE_FALSE_XOR_XNOR/"
name <- "slInZeros"
pdf(glue('SL_in_commonSea_pseudoAttractors_all_functions_noTRUE_FALSE_XOR_XNOR_{name}.pdf'))
res     <- list()
labels  <- list()
for (SLNUMBER in c(5,10,20))
{
    res[[length(res)+1]] <- read.table(glue('{path}{name}_{SLNUMBER}_augmAND.txt'))[,1]
    labels[[length(labels)+1]] <- glue('{SLNUMBER}_And')
    res[[length(res)+1]] <- read.table(glue('{path}{name}_{SLNUMBER}_augmOR.txt'))[,1]
    labels[[length(labels)+1]] <- glue('{SLNUMBER}_Or')
    res[[length(res)+1]] <- read.table(glue('{path}{name}_{SLNUMBER}_augmRND.txt'))[,1]
    labels[[length(labels)+1]] <- glue('{SLNUMBER}_Rnd')
}

boxplot(res, xaxt = "n", main="No. of selfloops in common sea with value 0", ylab="no. of selfloops")
axis(1, at = seq(1,length(res)), las = 2,labels = labels) # axis, ticks
dev.off()