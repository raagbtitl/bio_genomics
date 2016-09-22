
library(ggplot2); data1 <- read.table("common_slan_cd14.txt",header=T,sep="\t")

x<-suppressWarnings(as.numeric(as.character(data1[,2]))) #log2fc 
y<-suppressWarnings(-log10( as.numeric(as.character(data1[,4]))))  #pval 

##data1$color_flag <- ifelse(x > 1.5, 1, ifelse(x < -1.5,-1,0)) 
data1$threshold = as.factor(abs(x) > 1.5);xx<-na.omit(x);yy<-na.omit(y)
yy[!is.finite(yy)] <- 0 ;
g = ggplot(data=data1, aes(x=x, y=y,colour=threshold)) +  geom_point(size=1) + theme(legend.position = "none") + 
  geom_hline(yintercept = 1.30103,color="red",linetype="dashed")+
  xlim(c(min(xx), max(xx))) + ylim(c(0, max(yy)+100)) +   xlab("log2 fold change") + ylab("-log10 p-value");g



