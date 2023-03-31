library(glue)

type = "allFunctions" 
type = "noTRUE_FALSE_XOR_XNOR"
path=glue('n100k2p05_pseudoAttractors_{type}_29Marzo23')

pdf(glue('CORRELATIONS_{type}.pdf'))
par(mfrow=c(3,3), mai = c(0.58, 0.58, 0.17, 0.05))
for (SLNUMBER in c(5,10,20))
{
    for (SL_TYPE in c("augmAND", "augmOR", "augmRND"))
    {
        y.cs        <- read.table(glue('{path}/commonSeaSize_sl_{SLNUMBER}_{SL_TYPE}.txt'))[,1]
        y.cs_ones   <- read.table(glue('{path}/commonSeaOnes_sl_{SLNUMBER}_{SL_TYPE}.txt'))[,1]
        x.pa        <- read.table(glue('{path}/no_PseudoAttractors_{SLNUMBER}_{SL_TYPE}.txt'))[,1]
        
        #pdf(glue('CORRELATIONS_{type}_{SLNUMBER}_{SL_TYPE}_COMMON_SEA.pdf'))
        plot(x=x.pa, y=y.cs, 
                xlab="no. of (pseudo)attractors", 
                ylab= "Common sea size", 
                main=glue('{SLNUMBER}_{SL_TYPE}'),
                cex.main=0.9,
                pch=19, 
                col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5), 
                ylim=c(0,100))
        points(x=x.pa, y=y.cs_ones, 
                pch=19, 
                col=rgb(red = 0, green = 0, blue = 1, alpha = 0.2))
        legend("topright", 
            legend=c("tot. size", "sea = 1"), 
            col=c(rgb(red = 1, green = 0, blue = 0, alpha = 0.5), 
                rgb(red = 0, green = 0, blue = 1, alpha = 0.2)), 
            pch=19, 
            cex=0.8)
    }
}


dev.off()
