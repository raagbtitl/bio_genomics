#!/usr/bin/R
library(STRINGdb)
library( "org.Hs.eg.db" )

st<-read.csv("/home/rahul/ITLIVER/MACE/all_comparisons_MACE_ITliver.csv", header=T,sep="\t")

string_db <- STRINGdb$new( version="9_1", species=9606,score_threshold=0, input_directory="" )

ind <- with(st,(st$pvalue_TGFb_2h_vs_control_2h < 0.05) )

st1<-st[ind,] 

st2<-(st1[order(st1$pvalue_TGFb_2h_vs_control_2h),])

st2<-st2[1:50,]

mapped <- string_db$map( st2, "gene", removeUnmappedRows = TRUE )

hits <- mapped$STRING_id

mapped_pval05 <- string_db$add_diff_exp_color(mapped,logFcColStr="log2FoldChange_TGFb_2h_vs_control_2h" )

payload_id <- string_db$post_payload(mapped_pval05$STRING_id,colors=mapped_pval05$color )

##string_db$plot_network( hits, payload_id=payload_id )
enrichmentKEGG <- string_db$get_enrichment( hits, category = "KEGG", methodMT = "fdr", iea = TRUE )
##enrichmentGO <- string_db$get_enrichment( hits, category = "Process", methodMT = "fdr", iea = TRUE)
##png("stringdb_TGFb_2h_vs_control_2h")
write.table(enrichmentKEGG,"pathway_enriched_ppi_mace.csv",sep="\t")

pdf("stringdb_TGFb_2h_vs_control_2h_top50_up_down.pdf",width=6,height=4,paper='special') 
string_db$plot_network( hits, payload_id=payload_id )
dev.off()
