
files <- (list.files(pattern=".sbml"))
f <- function(bn) {return(tryCatch(loadSBML(bn),error=function(e) NULL))}
nets <- lapply(files, f)
nets[sapply(nets, is.null)] <- NULL #rimuove quelle con "problemi"

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
countSelfLoop <- function(bn) {sum(sapply(bn$genes,FUN=slFUNCTION, network=bn))}
countNodes <- function(bn) length(bn$genes)
fracSL_Nodes <- sapply(nets, FUN=function(bn) {countSelfLoop(bn)/countNodes(bn)})


"TBET" %in% reteprova$interactions[["TBET"]]$expression