
library(ReactomePA)
library(clusterProfiler)
library(knitr)
library( "org.Hs.eg.db" )
library( "reactome.db" )
library(DOSE)

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

##cleanTable <- table[ -which( is.na(table), arr.ind=TRUE )[,"row"], ]
##opts_chunk$set(tidy=TRUE,tidy.opts=list(keep.blank.line=FALSE, width.cutoff=50),out.truncate=80,out.lines=6,cache=TRUE,dev='pdf',include=TRUE,fig.width=6,fig.height=6,resolution=150)

##options(digits=3, width=80, prompt=" ", continue=" ")

dt<- read.csv("/home/rahul/ITLIVER/MACE/all_comparisons_MACE_ITliver.csv", header=T,sep="\t")

a<-as.vector(dt$gene)

dt$entrez <- convertIDs( a, "SYMBOL", "ENTREZID", org.Hs.eg.db )

ind <- with(dt,(dt$pvalue_LY_2h_vs_control_2h < 0.05) )

dt1<-dt[ind,] 

de1<- as.numeric(dt1$entrez[!is.na(dt1$entrez)])

x11 <- function( enzid ) {        ###perform ReactomePA
x1 <- enrichPathway(gene=enzid,organism = "human",pvalueCutoff=0.05, readable=T) ##hypergeometric test
write.table(summary(x1),file="pathway_mace_itl_LY_2h_vs_control_2h.csv",sep="\t",row.names=F)
barplot(x1, showCategory=10)
##cnetplot(x1, categorySize="pvalue",foldChange=dt$log2FoldChange_TGFb_2h_vs_control_2h)
}

x11(de1)  ###run ReactomePA


#### Execute reactome.db#####
reactomeTable <- AnnotationDbi::select( reactome.db, keys=dt1$entrez,keytype="ENTREZID", columns=c("ENTREZID","REACTOMEID","PATHNAME") )

incm <- do.call( rbind, with(reactomeTable, tapply( 
  ENTREZID, factor(REACTOMEID), function(x) dt1$entrez %in% x ) ))

colnames(incm) <- dt1$entrez
str(incm)

incm <- incm[ rowSums(incm) >= 5, ]

testCategory <- function( reactomeID ) {
  isMember <- incm[ reactomeID, ]
  
  data.frame( 
    reactomeID = reactomeID,
    numGenes = sum( isMember ),
    avgLFC = mean( dt1$log2FoldChange_LY_2h_vs_control_2h[isMember]),
    zValue = mean( dt1$log2FoldChange_LY_2h_vs_control_2h[isMember]) / sd( dt1$log2FoldChange_LY_2h_vs_control_2h[isMember]),
    strength = sum( dt1$log2FoldChange_LY_2h_vs_control_2h[isMember] ) / sqrt(sum(isMember)),
    pvalue = t.test( dt1$log2FoldChange_LY_2h_vs_control_2h[ isMember ])$p.value,
    reactomeName = reactomePATHID2NAME[[reactomeID]] ) 
}

reactomeResult <- do.call( rbind, lapply( rownames(incm), testCategory ) )

reactomeResult$padjust <- p.adjust( reactomeResult$pvalue, "BH" )

write.csv(reactomeTable,file="pathway_genes_LY_2h_vs_control_2h.csv", row.names=F)

write.csv(reactomeResult,file="pathway_result_score_LY_2h_vs_control_2h.csv", row.names=F)


