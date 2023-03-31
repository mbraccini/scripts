library(glue)

type = "allFunctions" 
type = "noTRUE_FALSE_XOR_XNOR"
path=glue('n100k2p05_pseudoAttractors_{type}_29Marzo23')


for (SLNUMBER in c(5,10,20))
{
    for (SL_TYPE in c("augmAND", "augmOR", "augmRND"))
    {
        y.cs <- read.table(glue('{path}/commonSeaSize_sl_{SLNUMBER}_{SL_TYPE}.txt'))[,1]
        x.pa <- read.table(glue('{path}/no_PseudoAttractors_{SLNUMBER}_{SL_TYPE}.txt'))[,1]
        print(x)
        print(length(x))
        print(y)
        print(length(y))
        pdf(glue('CORRELATIONS_{type}_{SLNUMBER}_{SL_TYPE}.pdf'))
        plot(x=x.pa, y=y.cs, xlab="no. of (pseudo)attractors", ylab= "Common sea size", main=glue('{type}_{SLNUMBER}_{SL_TYPE}'))
        dev.off()
    }
}


