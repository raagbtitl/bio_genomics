require(ConsensusClusterPlus)
require(pamr)
rem <- function(x){   x <- t(apply(x,1,as.numeric))
                      r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
                      remove <- which(r > dim(x)[2]*0.5)
                      return(remove) }

d<-read.csv("LIHC.rnaseqv2_RSEM_genes_normalized_data.csv",header=T,dec=".",sep=",",row.names=1);d<-as.matrix(d)
remove <- rem(d)
d <- d[-remove,]
mads=apply(d,1,mad); 
d=d[order(mads,decreasing=T)[1:3000],]
d = sweep(d,1, apply(d,1,median,na.rm=T))

results<-ConsensusClusterPlus(d,maxK=8,reps=1000,pItem=0.8,pFeature=1,title="lihc_pearson_hc",distance="pearson",clusterAlg="hc",plot="png")
table(results[[5]]$consensusClass)

findDiffGenes <- function(X, groups, top=500){
  f <- factor(groups)
  topGenes <- lapply(levels(f), function(lvl){
    gf <- factor(lvl == f)
    adjp <- p.adjust(apply(X, 1, function(x){ wilcox.test(x ~ gf)$p.value}), method="BH")
    # just take top genes, rather than filter on adjusted pvalue
    return (rownames(X)[order(adjp)[1:top]])
  })
  unique(do.call("c", topGenes))
}


diffGenes <- findDiffGenes(d, results[[5]]$consensusClass)
ds <- d[rownames(d) %in% diffGenes,]

par(mfrow=c(2,2))
pc1 <- svd(d - rowMeans(d))
plot((pc1$d^2 / sum(pc1$d^2))[1:100],main="Eig distribution",ylab="% var")
plot(pc1$v[,1], pc1$v[,2], col=factor(results[[5]]$consensusClass), xlab="PC1",ylab="PC2",
     main="PC Plot (Pre Gene Selection)",pch=19,cex=.8)

pc2 <- svd(ds - rowMeans(ds))
plot((pc2$d^2 / sum(pc2$d^2))[1:100],main="Eig distribution",ylab="% var")
plot(pc2$v[,1], pc2$v[,2], col=factor(results[[5]]$consensusClass), xlab="PC1",ylab="PC2",
     main="PC Plot (Post Gene Selection)",pch=19,cex=.8)
# build classifier


labels <- as.factor(as.numeric(results[[5]]$consensusClass))

sample=colnames(ds)
geneid<-rownames(ds)

train <- list(x=ds, y=labels ,genenames = geneid, geneid = geneid ,sampleid = sample)

fit <- pamr.train(data=train)

cvfit <- pamr.cv(fit, train,nfold=10)

opt.threshold <- cvfit$threshold[which(cvfit$error == min(cvfit$error))]
if(length(opt.threshold) > 1){
  opt.threshold <- opt.threshold[length(opt.threshold)]
}

pamr.plotcv(cvfit)

Delta = 4 ;

pamr.confusion(cvfit, Delta)

pamr.plotcen(fit, train, Delta)

pamr.plotcvprob(fit, train, Delta)

pamr.geneplot(fit, train, Delta)

pamr.listgenes(fit, train, Delta, genenames = TRUE)

thresh <- fit$threshold[18]

Xmean <- rowMeans(ds)
Xsd <- apply(ds, 1, sd)

predList <- lapply(ds, function(eset){
  X <-ds
  idxs <- match(rownames(ds), rownames(X))
  Xtest <- X[idxs,]
  pred <- pamr.predict(fit, Xtest,  opt.threshold, type="posterior")
  return (pred)
})

write.table(pred, file="pred.tsv",quote=FALSE,row.names=FALSE)
