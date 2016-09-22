require(gage)
library( "org.Hs.eg.db" )

kg.hsa <- kegg.gsets( "hsa" )
kegg.gs2 <- kg.hsa$kg.sets[ kg.hsa$sigmet.idx ]

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

dt<-na.omit(dt)

ind <- with(dt,{(dt$log2FoldChange_TGFb_2h_vs_control_2h < -1.5) | (dt$log2FoldChange_TGFb_2h_vs_control_2h > 1.5)})

dt<-dt[ind,] 

rownames(dt)<-dt$entrez

fc<-dt$log2FoldChange_TGFb_2h_vs_control_2h

names(fc)=rownames(dt)

res.kegg <- gage( fc, gsets= kegg.gs2,ref= NULL,samp= NULL)

write.table(res.kegg,file="kegg_pathwayTGFb_2h_vs_control_2h.csv", row.names = T,sep='\t')

sel <- res.kegg$greater[, "p.val"] < 0.05 & !is.na(res.kegg$greater[,"p.val"])

path.ids <- rownames(res.kegg$greater)[sel]
path.ids

##path.ids2 <- substr(path.ids, 1, 8)

##pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = names(fc), pathway.id = pid, species = "hsa"))



