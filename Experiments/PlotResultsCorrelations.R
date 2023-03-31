library(glue)

type = "allFunctions" 
type = "noTRUE_FALSE_XOR_XNOR"
path=glue('n100k2p05_pseudoAttractors_{type}_29Marzo23')


for (SLNUMBER in c(5,10,20))
{
    for (SL_TYPE in c("augmAND", "augmOR", "augmRND"))
    {
        y.cs        <- read.table(glue('{path}/commonSeaSize_sl_{SLNUMBER}_{SL_TYPE}.txt'))[,1]
        y.cs_ones   <- read.table(glue('{path}/commonSeaOnes_sl_{SLNUMBER}_{SL_TYPE}.txt'))[,1]
        x.pa        <- read.table(glue('{path}/no_PseudoAttractors_{SLNUMBER}_{SL_TYPE}.txt'))[,1]
        print(x)
        print(length(x))
        print(y)
        print(length(y))
        pdf(glue('CORRELATIONS_{type}_{SLNUMBER}_{SL_TYPE}_COMMON_SEA.pdf'))
        plot(x=x.pa, y=y.cs, xlab="no. of (pseudo)attractors", ylab= "Common sea size", main=glue('{type}_{SLNUMBER}_{SL_TYPE}'),pch=19, col=rgb(red = 1, green = 0, blue = 0, alpha = 0.4), ylim=c(0,100))
        points(x=x.pa, y=y.cs_ones, pch=19, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.1))
        legend("topright", 
            legend=c("Common sea size", "Common sea with value = 1"), 
            col=c(rgb(red = 1, green = 0, blue = 0, alpha = 0.4), rgb(red = 0, green = 0, blue = 1, alpha = 0.1)), 
            pch=19, 
            cex=0.8)
        dev.off()
    }
}


