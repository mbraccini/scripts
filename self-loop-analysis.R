
library(BoolNet)
library(glue)
loadSBMLnetworks <- function(){
	files <- (list.files(pattern=".sbml"))
	f <- function(bn) {return(tryCatch(loadSBML(bn),error=function(e) NULL))}
	nets <- lapply(files, f)
	nets[sapply(nets, is.null)] <- NULL #rimuove quelle con "problemi"
	return (nets)
}
loadNetworks <- function(myPath, myPattern){
	files <- list.files(path=myPath, pattern = myPattern)
	files <- paste0(myPath,files)
	f <- function(bn) {return(tryCatch(BoolNet::loadNetwork(bn),error=function(e) NULL))}
	nets <- lapply(files, f)
	nets[sapply(nets, is.null)] <- NULL #rimuove quelle con "problemi"
	return(nets)
}

#countSelfLoop <- function(bn) {sum(sapply(bn$genes, FUN=function(geneName) grepl(geneName, bn$interactions[[geneName]]$expression)))}
#countSelfLoop <- function(bn) {sum(sapply(bn$genes,FUN=function(geneName) {grepl(geneName, bn$interactions[[geneName]]$expression) && !(geneName %in% bn$interactions[[geneName]]$expression)}))}
slFUNCTION <- function(geneName, network) {

	expr <- network$interactions[[geneName]]$expression
	exprREMOVED <- gsub("[[:blank:]()!]", "", expr, fixed = FALSE)
	exprSPLITTED <- unlist(strsplit(exprREMOVED,"[|&]"))

	if (length(exprSPLITTED) == 1){
		return(FALSE)
	} else {
		return (geneName %in% exprSPLITTED)
	}

}

nets <- loadNetworks("./Experiments/testing/", "_1$")

countSelfLoop <- function(bn) {sum(sapply(bn$genes,FUN=slFUNCTION, network=bn))}
countNodes <- function(bn) length(bn$genes)
fracSL_Nodes <- sapply(nets, FUN=function(bn) {countSelfLoop(bn)/countNodes(bn)})

#"TBET" %in% reteprova$interactions[["TBET"]]$expression


folder="conteggioAutoanelli_boolnet"
dir.create(folder)

samples <- 100
stats <- list()
startNodes <- 10
endNodes <- 100
delta <- 10
for (N in seq(startNodes,endNodes,delta)){

	for (i in seq(1,samples)){
		res <- vector("integer", length = samples)

		bn <- generateRandomNKNetwork(n=N, 
										k=2,
										functionGeneration = "biased", 
										zeroBias = 0.5,
										readableFunctions=TRUE,
										noIrrelevantGenes=FALSE)
		#print(bn$interactions)
		BoolNet::saveNetwork(bn, glue('{folder}/{N}_{i}'))
		#rete <-BoolNet::loadNetwork(glue('{folder}/{N}_{i}'))
		res[i] <- countSelfLoop(bn)/countNodes(bn)
	}	
	stats[[length(stats) +1 ]] <- res
}
print(stats)
pdf("BoolNetFractionSelfloops.pdf")
boxplot(stats, xaxt = "n", main="Selfloops in RBN (BoolNet)", ylab="Fraction of selfloops")
axis(1, at = seq(1,((endNodes-startNodes)/delta)+1), labels = seq(startNodes, endNodes,by=delta))
dev.off()
