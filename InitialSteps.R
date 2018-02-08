rna <- rna <- read.table("KIRP.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",header=T,sep="\t",check.names = FALSE);
rna<-rna[!duplicated(rna[,1]),]; 
rownames(rna)<-rna[,1];rna<-rna[,-1]
rna<-as.matrix(rna)
is.numeric(rna)
table(substr(colnames(rna),14,14))

rna_cp<-rna

# get the index of the normal/control samples
n_index <- which(substr(colnames(rna),14,14) == "1")
t_index <- which(substr(colnames(rna),14,14) == "0")

# first get rid of genes whose expression is == 0 in more than 50% of the samples:
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i ==0)))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}
remove <- rem(rna)
rna<-rna[-remove,]
dim(rna)
is.numeric(rna)

#  genes whose expression is gt 5 in more than 80% of the samples:
kept <- function(rna){
  rna<- as.matrix(rna)
  rna<- t(apply(rna,1,as.numeric))
  r <- as.numeric(apply(rna,1,function(i) sum(i >5)))
  keep <- which(r > dim(rna)[2]*0.8)
  return(keep)
}
keep <- kept(rna)
rna<-rna[keep,]
dim(rna)
is.numeric(rna)

selTopVarGenesByCV <- function(rna){                        ##Coeff of Var based filtering
  ffun=filterfun(cv(a = 0.8,b=10))
  ##pOverA(p = 0.1, A = 100), 
  filt=genefilter(rna,ffun); sdat<-rna[filt,]
  return(sdat)
}

selTopVarGenesByMAD <- function(rna){                        ##MAD based filtering
  # select genes with highest mad (quantile > 60%)
  madrna <- apply(rna, 1, mad)
  rna <- rna[madrna > quantile(madrna, .8),]
  return(rna)
}
selTopVarGenesByCVq3 <- function(rna){
  mean <- abs(rowMeans(rna))  # mean of normal
  sd <- apply(rna,1,sd)  # SD of normal
  cv<-sd/mean
  q<-ceiling(length(cv)*0.8)
  sortcv<-sort(cv)
  cvfun<-which(cv>sortcv[q])
  sdat<-rna[cvfun,]
  return(sdat)
}

rna <- selTopVarGenesByCVq3(rna)
dim(rna)
colnames(rna) <- gsub("\\.","-",substr(colnames(rna),1,12))

distances <- c("euclidean", "pearson", "spearman")

calculate_distance <- function(data, method) {
  return(if (method == "spearman") {
    as.matrix(1 - cor(data, method = "spearman"))
  } else if (method == "pearson") {
    as.matrix(1 - cor(data, method = "pearson"))
  } else {
    as.matrix(dist(t(data),method="euclidean"))
  })
}

dists<-calculate_distance(rna,distances[1])

transformation <- function(dists, method) {
  if (method == "pca") {
    t <- prcomp(dists, center = TRUE, scale. = TRUE)
    return(t$rotation)
  } else if (method == "laplacian") {
    L <- norm_laplacian(dists)
    l <- eigen(L)
    # sort eigenvectors by their eigenvalues
    return(l$vectors[, order(l$values)])
  }
}

transformed_rna<-transformation(dists,"pca")

# Check for the optimal number of clusters given the data using Elbow method

mydata <- transformed_rna
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares",
     main="Assessing the Optimal Number of Clusters with the Elbow Method",
     pch=20, cex=2)


# Perform K-Means with 2 clusters
set.seed(7)
i<-     ##based on above observation
km1 = kmeans(transformed_rna, i, nstart=100)

# Plot results
plot(dat, col =(km1$cluster +1) , main="K-Means result with 2 clusters", pch=20, cex=2)

fileWithClsNum<-read.csv("patients with cluster number.csv", sep=",", row.names = 1,stringsAsFactors = F)

ks<- c(1:5)

for (i in ks){
subClusterPat<-which(fileWithClsNum[,2]==i)
}

 ############  for each sub cluster run differnetial expression analysis between live and dead patients ############
status<-read.table("A_vs_D.csv",sep="\t",row.names=1)
df<-sapply(status["patient.bcr_patient_barcode",], toupper)
a_index <- which(status[2,]== "alive" & status[3,]== "tumor free")
a_patient<-df[a_index]
d_index <- which( status[2,]== "dead" & status[3,]== "with tumor")
d_patient<-df[d_index]

colnames(rna_cp) <- gsub("\\.","-",substr(colnames(rna_cp),1,12))
nodupind<-which(table(colnames(rna_cp))<=1)
rna_vm<-rna_vm[,nodupind]
rnanew<-cbind(rna_cp[,colnames(rna_cp) %in% a_patient ],rna_cp[,colnames(rna_cp) %in% d_patient])

train <- t(rnanew)
labs <- factor(rownames(train))
rownames(train) <- NULL
model <- tryCatch(svm(train, labs, kernel = kern), error = function(cond) return(NA))
pred <- predict(model, t(study))

pca <- as.data.frame(prcomp(t(na.omit(rnanew)))$x)
plot(PC2~PC1, data=pca, col=cond, pch=19, main="PC Plot")

##### Run DE Analysis between living and dead patients using limma package #####

library(limma)

cond <- factor(ifelse(seq(1,dim(rnanew)[2],1) %in% which(colnames(rnanew) %in% a_patient), 1,  0))
d <- model.matrix(~1+cond)

fit <- eBayes(lmFit(rnanew,d)) 
limma::plotMA(fit); 

res=topTable(fit,coef=2,n=nrow(fit$t))  

res[,7]<-exp(res[,6])
res[,8]<-res[,7]/(1+res[,7])
res[,9] <- ifelse((res$P.Val < 0.01 & abs(res$logFC) > 1.5), "red", "black") 
size <- ifelse((res$P.Val < 0.01 & abs(res$logFC) > 1.5), 4, 2)

##Construct the plot object
require(ggplot2)

g <- ggplot(data=res, aes(x=res[,1], y=-log10(res$P.Val))) +
  geom_point(size=size, colour=res[,9]) +
  xlim(c(-3, 3)) + ylim(c(0,8)) +
  xlab("log2 fold change") + ylab("-log10  p-value") +
  guides(colour = guide_legend(override.aes = list(shape=16)))

g

#### sc3######
get_de_genes <- function(dataset, labels) {
  tmp <- apply(dataset, 1, kruskal.test, g = factor(labels))
  ps <- unlist(lapply(tmp, "[[", "p.value"))
  ps <- p.adjust(ps)
  return(ps)
}

support_vector_machines <- function(train, study, kern) {
  train <- t(train)
  labs <- factor(rownames(train))
  rownames(train) <- NULL
  model <- tryCatch(svm(train, labs, kernel = kern), error = function(cond) return(NA))
  pred <- predict(model, t(study))
  return(pred = pred)
}

reindex_clusters <- function(hc, k) {
  clusts <- stats::cutree(hc, k)
  labels <- names(clusts)
  names(clusts) <- 1:length(clusts)
  ordering <- clusts[hc$order]
  new.index <- NULL
  j <- 1
  for (i in unique(ordering)) {
    tmp <- rep(j, length(ordering[ordering == i]))
    names(tmp) <- names(ordering[ordering == i])
    new.index <- c(new.index, tmp)
    j <- j + 1
  }
  clusts <- new.index
  clusts <- clusts[order(as.numeric(names(clusts)))]
  names(clusts) <- labels
  return(clusts)
}
#### Estimate the optimal k for k-means clustering based on number of 
#### significant eigenvalues according to the Tracy-Widom test
#### Can use pkg -- https://github.com/patperry/r-rmtstat or import function from SC3

estkTW <- function(dataset) {
  
  p <- ncol(dataset)
  n <- nrow(dataset)
  
  # compute Tracy-Widom bound
  x <- scale(dataset)
  muTW <- (sqrt(n - 1) + sqrt(p))^2
  sigmaTW <- (sqrt(n - 1) + sqrt(p)) * (1/sqrt(n - 1) + 1/sqrt(p))^(1/3)
  sigmaHatNaive <- tmult(x)  # x left-multiplied by its transpose
  bd <- 3.273 * sigmaTW + muTW  # 3.2730 is the p=0.001 percentile point for the Tracy-Widom distribution
  
  # compute eigenvalues and return the amount which falls above the bound
  evals <- eigen(sigmaHatNaive, symmetric = TRUE, only.values = TRUE)$values
  k <- 0
  for (i in 1:length(evals)) {
    if (evals[i] > bd) {
      k <- k + 1
    }
  }
  return(k)
}




