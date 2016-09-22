library(topGO)

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

a<-as.vector(dt$gene)

dt$entrez <- convertIDs( a, "SYMBOL", "ENTREZID", org.Hs.eg.db )

go2entrez <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "Entrez")

allGenes <- unique(unlist(go2entrez ))

##fdr <- p.adjust(dt$pvalue_TGFb_2h_vs_control_2h, method = "BH");names(fdr) <- dt$entrez
##f <- function(q) {  return (q < 0.01)}
##topgo <- new("topGOdata",ontology = "CC", allGenes = fdr, geneSel = f,annot = annFUN.org,mapping = "org.Hs.eg.db",ID = "Entrez")

geneList <- factor(as.integer(allGenes %in% dt$entrez))

names(geneList)<-allGenes

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel = function(p) p < 0.01, description = "Test", annot = annFUN.org, mapping = "org.Hs.eg.db",ID = "Entrez")

resultFisher <- runTest(GOdata, algorithm='classic', statistic='fisher')

sel.terms <- usedGO(GOdata)

allRes <- GenTable(GOdata, classic=resultFisher, topNodes=length(sel.terms), numChar=9999999)

ann.genes <- genesInTerm(GOdata, sel.terms)

ag<-as.data.frame(stack(ann.genes));agf<-ag[,c(2,1)]

agf2<-as.data.frame(cbind(as.character.factor(agf$ind),agf$values));colnames(agf2)<-c("GO.ID", "gene")

mdata<-aggregate(agf2$gene ~ agf2$GO, data=agf2, FUN=paste, collapse=",")

mdata1<-as.data.frame(cbind(mdata["agf2$GO"],mdata["agf2$gene"]))

colnames(mdata1)<-c("GO.ID", "gene")

mrge.data<-merge(allRes,mdata1,by="GO.ID")


##write.table(agf2,quote = FALSE,col.names = F,row.names = F,file=paste("BP","_GO_transcript_annots.csv", sep=""))
##write.table(mdata,quote = FALSE,col.names = F,row.names = F,file=paste("BP","_GO_transcript_annots.slimmed.csv", sep=""))
##write.table(allRes, file="topGo_results.1.csv", sep="\t", row.names=FALSE)
write.table(mrge.data,quote = FALSE,col.names = F,row.names = F,file=paste("BP","_GO_transcript_annots.slimmed.enzid.csv", sep=""))