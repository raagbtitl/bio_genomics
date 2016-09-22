
library( "org.Hs.eg.db" )

library( "reactome.db" )

dt<-read.csv("LPA Hippo targets_RA-1.csv",header=T)
  ##read.table("gene_path.txt",header=T)


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

a<-as.vector(dt$gene)

dt$entrez <- convertIDs( a, "SYMBOL", "ENTREZID", org.Hs.eg.db )

reactomeTable <- AnnotationDbi::select( reactome.db, keys=dt$entrez, 
                                         keytype="ENTREZID", columns=c("ENTREZID","REACTOMEID","PATHNAME") )

reactomeTable$Ens <- convertIDs( reactomeTable$ENTREZID, "ENTREZID", "ENSEMBL",org.Hs.eg.db )

reactomeTable$gene <- convertIDs( reactomeTable$Ens, "ENSEMBL", "SYMBOL",org.Hs.eg.db )

incm <- do.call( rbind, with(reactomeTable, tapply( 
     ENTREZID, factor(REACTOMEID), function(x) dt$entrez %in% x ) ))

colnames(incm) <- dt$entrez
str(incm)

incm <- incm[ rowSums(incm) >= 5, ]


testCategory <- function( reactomeID, dt$log2FoldChange_TGFb_2h_vs_control_2h ) {
     isMember <- incm[ reactomeID, ]
     log2FoldChange=dt$log2FoldChange_TGFb_2h_vs_control_2h
     data.frame( 
         reactomeID = reactomeID,
         numGenes = sum( isMember ),
         avgLFC = mean( log2FoldChange[isMember] ),
         zValue      = mean( log2FoldChange[isMember] ) / sd( log2FoldChange[isMember] ),
         strength = sum( log2FoldChange[isMember] ) / sqrt(sum(isMember)),
         pvalue = t.test( log2FoldChange[ isMember ] )$p.value,
         reactomeName = reactomePATHID2NAME[[reactomeID]] ) }


reactomeResult <- do.call( rbind, lapply( rownames(incm), testCategory ) )

reactomeResult$padjust <- p.adjust( reactomeResult$pvalue, "BH" )

##reactomeResultSignif <- reactomeResult[ reactomeResult$padjust < 0.05, ]

##reactomeResultSignif[ order(reactomeResultSignif$strength), ]

write.csv(reactomeResult,file="pathway_result.csv", row.names=F)

write.csv(reactomeTable,file="pathway_genes.csv", row.names=F)

##e2s = toTable(org.Hs.egSYMBOL)
##e2s$symbol[e2s$gene_id %in% reactomeTable$ENTREZID]