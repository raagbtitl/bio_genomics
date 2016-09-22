library(TCGA2STAT)
library(CCA)
rnaseq2 <- getTCGA(disease="LIHC", data.type="RNASeq2")
methyl <- getTCGA(disease="LIHC", data.type="Methylation", type="450K",clinical=FALSE)

# Use only the genes with most perturbed expression and mehylation
met.var <- apply(methyl[,-c(1,2,3)], 1, var)
met.data <- subset(methyl[,-c(1,2,3)], met.var >= quantile(met.var, 0.99, na.rm=T)
                                    & !is.na(met.var))
                 
rnaseq2.var <- apply(log10(1+rnaseq2$dat), 1, var)
rnaseq.data <- subset(log10(1+rnaseq2$dat), rnaseq2.var >=quantile(rnaseq2.var, 0.99, na.rm=T) & !is.na(rnaseq2.var))
                 
# The package utility function merges the two data of same patients
met.rnaseq2 <- GeneMerge(dat1 = rnaseq.data, dat2= met.data)
                 
 # Carried out CCA on X and Y 
estl <- estim.regul(met.rnaseq2$X, met.rnaseq2$Y)
                 
lusc.cc <- rcc(met.rnaseq2$X, met.rnaseq2$Y, estl$lambda1, estl$lambda1)
plt.cc(lusc.cc, d1=1, d2=2, type="b", var.label=TRUE)