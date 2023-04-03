library(glue)


#type = "allFunctions" 
type = "noTRUE_FALSE_XOR_XNOR"
path_100=glue('n100k2p05_pseudoAttractors_{type}_31Marzo23')
path_300=glue('n300k2p05__{type}_noPseudoAttsAsSizeIncreases_3Aprile')
path_500=glue('n500k2p05__{type}_noPseudoAttsAsSizeIncreases_3Aprile')
path_700=glue('n700k2p05__{type}_noPseudoAttsAsSizeIncreases_3Aprile')



#expString   <- "no_attractors_"
expString   <- "no_PseudoAttractors_"
pdf(glue('{type}_{expString}AsSizeIncreases.pdf'))

res         <- list()
labels      <- list()
SLNUMBER    <- 20

#100
res[[length(res)+1]] <- read.table(glue('{path_100}/{expString}{SLNUMBER}_augmAND.txt'))[,1]
labels[[length(labels)+1]] <- glue('100_And')
res[[length(res)+1]] <- read.table(glue('{path_100}/{expString}{SLNUMBER}_augmRND.txt'))[,1]
labels[[length(labels)+1]] <- glue('100_Rnd')
#300
res[[length(res)+1]] <- read.table(glue('{path_300}/{expString}{SLNUMBER}_augmAND.txt'))[,1]
labels[[length(labels)+1]] <- glue('300_And')
res[[length(res)+1]] <- read.table(glue('{path_300}/{expString}{SLNUMBER}_augmRND.txt'))[,1]
labels[[length(labels)+1]] <- glue('300_Rnd')
#500
res[[length(res)+1]] <- read.table(glue('{path_500}/{expString}{SLNUMBER}_augmAND.txt'))[,1]
labels[[length(labels)+1]] <- glue('500_And')
res[[length(res)+1]] <- read.table(glue('{path_500}/{expString}{SLNUMBER}_augmRND.txt'))[,1]
labels[[length(labels)+1]] <- glue('500_Rnd')
#700
res[[length(res)+1]] <- read.table(glue('{path_700}/{expString}{SLNUMBER}_augmAND.txt'))[,1]
labels[[length(labels)+1]] <- glue('700_And')
res[[length(res)+1]] <- read.table(glue('{path_700}/{expString}{SLNUMBER}_augmRND.txt'))[,1]
labels[[length(labels)+1]] <- glue('700_Rnd')

boxplot(res, 
        xaxt = "n", 
        log='y', 
        ylab=glue('{expString}'), 
        xlab="nodes", 
        main="RBN with 20 selfloops augm{AND,RND}",
        col=c('lightblue','orange'),
        boxwex=.7)

axis(1, 
    at = seq(1,length(res)), 
    las = 2,
    labels=labels) # axis, ticks

dev.off()