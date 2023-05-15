library(glue)
library("gplots")
library(dplyr)
#type = "allFunctions" 
type = "noTRUE_FALSE_XOR_XNOR"
path=glue('n100k2p05_pseudoAttractors_{type}_26Aprile23')

SLTYPE <- "_augmAND"#"_augmOR"

for (SLNUMBER in c(1,2,3,4,5,10,20))
{
    NET <- 1
    specific_path <- glue('{path}/sl{SLNUMBER}{SLTYPE}/atts/bn_{NET}_pseudoAtts.csv')
    df <- read.csv(specific_path, header=TRUE)
    pdf(glue('prova_bn_{NET}_sl{SLNUMBER}.pdf'))
    #df <- data.frame(c(0,0,1),c(0,0,1),c(0,1,1))
    df <- df %>% arrange(across(everything())) 
    heatmap.2(as.matrix(df),
            dendrogram='none', 
            Rowv=FALSE, 
            Colv=FALSE,
            trace='none',
            scale="none", 
            main=glue('bn_{NET} with {SLNUMBER} selfloops'), 
            col=grey(seq(0, 1, length = 256)),
            key=FALSE, 
            density.info="none"
    )

    #image(t(df),col=grey(seq(0, 1, length = 256)),xaxt='n',yaxt='n')
    #heatmap(as.matrix(df))
    dev.off()
    
}
    