setwd("~/Desktop")
data<-read.csv("1_CT_CPM_promoter_adj.csv",header=T,sep="\t",check.names=F,row.names=1)
plotMA<-function(res)
{
x1<-sqrt(res$baseMeanA*res$baseMeanB);x<-log2(x1)
##res$baseMean
plot(x, res$log2FoldChange,col = ifelse(res$pval < 0.05, "red", "black"), ylim=c(-5,5), main = "",pch=16, cex=0.7, log = "x", frame.plot=FALSE,ylab="Log2 fold change",xlab="Log2 normalized mean expression")
}

plotMA(data)
