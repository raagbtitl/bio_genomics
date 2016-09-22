install.packages("ggplot2")

gene_list <- read.table("/Users/Javi/Desktop/gene_list.csv", header=T, sep=",")

require(ggplot2)
##Highlight genes that have an absolute fold change > 2 and a p-value < 0.05
gene_list$threshold = as.factor(abs(gene_list$logFC) > 2 & gene_list$P.Value < 0.05)

##Construct the plot object
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value), colour=my_palette)) +
  geom_point(alpha=0.4, size=5) +
  theme(legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
g
gene_list$color_flag <- ifelse(gene_list$logFC > 1.3, 1, ifelse(gene_list$logFC < -1.3, -1, 0))
##gene_list$color_flag <- ifelse(gene_list$logFC > 1.3, 1, ifelse(gene_list$logFC < -1.3, -1, 0))