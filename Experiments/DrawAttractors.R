library(glue)
library("gplots")
library(dplyr)

myplot <- function(df, filename, netnumber, slnumber, sltype){
    pdf(filename)
    #df <- data.frame(c(0,0,1),c(0,0,1),c(0,1,1))
    df <- df %>% arrange(across(everything())) 
    heatmap.2(as.matrix(df),
            dendrogram='none', 
            Rowv=FALSE, 
            Colv=FALSE,
            trace='none',
            scale="none", 
            main=glue('bn_{netnumber} with {slnumber} selfloops {sltype}_'), 
            col=grey(seq(0, 1, length = 256)),
            key=FALSE, 
            density.info="none"
    )

    #image(t(df),col=grey(seq(0, 1, length = 256)),xaxt='n',yaxt='n')
    #heatmap(as.matrix(df))
    dev.off()
}

#type = "allFunctions" 
type = "noTRUE_FALSE_XOR_XNOR"
path=glue('n100k2p05_pseudoAttractors_{type}_26Aprile23')


NET <- 1

SLTYPE <- "_augmAND"#"_augmOR"
folder.path = glue("{path}/atts_plots/{gsub('^.', '', SLTYPE)}/")
dir.create(folder.path)

path_0sl <- glue('{path}/sl0/atts/bn_{NET}_pseudoAtts.csv')
df <- read.csv(path_0sl, header=TRUE)
myplot(df,glue('{folder.path}bn_{NET}_sl0{SLTYPE}.pdf'), NET, 0,"")

for (SLNUMBER in c(1,2,3,4,5,10,20))
{
    specific_path <- glue('{path}/sl{SLNUMBER}{SLTYPE}/atts/bn_{NET}_pseudoAtts.csv')
    df <- read.csv(specific_path, header=TRUE)
    myplot(df,glue('{folder.path}bn_{NET}_sl{SLNUMBER}{SLTYPE}.pdf'), NET, SLNUMBER, SLTYPE)
}
    