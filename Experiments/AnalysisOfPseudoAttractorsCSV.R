library(glue)

#type = "allFunctions" 
type = "noTRUE_FALSE_XOR_XNOR"
path=glue('n100k2p05_pseudoAttractors_{type}_26Aprile23')

res <- list()
#names(a)[1] <- c("biill")
for (SLNUMBER in c(1,2,3,4,5,10,20))
{
    for (SLTYPE in c("_augmAND","_augmOR")){
        numberOfAttractors <- c()
        for (NET in seq(1,100)){
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

        }
        res[[length(res) + 1]] <- numberOfAttractors
        names(res)[length(res)] <- glue('{SLNUMBER}{SLTYPE}')
    }
}

print(res)
print(names(res))

pdf(glue('filtered_noPseudoAttrs_{type}.pdf'))
par( mar = c(6.5, 4, 2, 2))

boxplot(res, xaxt = "n",  ylab="no. of (pseudo)attractors")
axis(1, at = seq(1,length(res)), las = 2,labels = names(res)) # axis, ticks
dev.off()

#pdf(glue('DISTANCES_{type}_sl_0.pdf'))
#dev.off()
