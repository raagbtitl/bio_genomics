library(GOstats)
library( "org.Hs.eg.db" )
convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select( 
    db, keys=ids, keytype=fromKey, cols=c(fromKey,toKey) ) )
  if( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

dt<- read.csv("/home/rahul/ITLIVER/MACE/all_comparisons_MACE_ITliver.csv", header=T,sep="\t")
a<-as.vector(dt$gene);dt$entrez <- convertIDs( a, "SYMBOL", "ENTREZID", org.Hs.eg.db )

selected <- unique(dt$entrez)

param <- new("GOHyperGParams", geneIds=selected,
             #universe=universe,
             annotation="org.Hs.eg.db", ontology="BP",pvalueCutoff=0.05,
             conditional=FALSE, testDirection="over")

hyp <- hyperGTest(param)
sumTable <- summary(hyp)
# subset the output table to get the columns of interest 
# (GO ID, GO description, p-value)
out <- subset(sumTable, select=c(1, 7, 2))
# retrieve input genes associated with each GO identifier
# use the org.Hs.eg data mapping to get GO terms for each ID
goMaps <- lapply(out$GOBPID, function(x) unlist(mget(x, org.Hs.egGO2ALLEGS)))
# subset the selected genes based on those in the mappings
goSelected <- lapply(goMaps, function(x) selected[selected %in% x])
# join together with a semicolon to make up the last column
out$inGenes <- unlist(lapply(goSelected, function(x) paste(x, collapse=";")))
# write the final data table as a tab separated file 
write.table(out, file="go_results.csv", sep="\t", row.names=FALSE)