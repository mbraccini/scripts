library(glue)

type = "allFunctions" 
#type = "noTRUE_FALSE_XOR_XNOR"
path=glue('n100k2p05_pseudoAttractors_{type}_31Marzo23')

pdf(glue('DISTANCES_{type}.pdf'))
par(mfrow=c(3,3), mai = c(0.53, 0.55, 0.14, 0))

for (SLNUMBER in c(5,10,20))
{
    for (SL_TYPE in c("augmAND", "augmOR", "augmRND"))
    {
        file=glue('{path}/distances_sl_{SLNUMBER}_{SL_TYPE}.csv')
        df<- read.table(file,sep=",",fill=TRUE)
        ls <- unlist(df)
        distances <- ls[!is.na(ls)]
        #hist(distances, breaks=0:(max(distances)+1), right=FALSE, freq=FALSE)
        tmp <- hist(distances, 
                    breaks=0:(max(distances)+1), 
                    xaxt="n", 
                    right=FALSE, 
                    freq=TRUE,
                    cex.lab=1, 
                    cex.axis=1.1,
                    cex.main=0.9,
                    ylab="Absolute Frequency",
                    xlab="Hamming distance",
                    main=glue('{SLNUMBER}_{SL_TYPE}'))
        #tmp$counts=tmp$counts/sum(tmp$counts)
        axis(1, at=tmp$mids, labels=0:max(distances),cex.axis=1.1)

    }
}


dev.off()
