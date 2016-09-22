

library(parallel)

run.edger <- function(counts, labels, normmethod="TMM",reorder=TRUE) {
  require(edgeR, quietly=TRUE)
  d <- DGEList(counts=counts, group=labels)
  if (normmethod != "none")
    d <- calcNormFactors(d, method=normmethod)
  d <- estimateCommonDisp(d)
  d <- estimateTagwiseDisp(d)
  et <- exactTest(d)
  p <- et$table$PValue
  adj.p <- p.adjust(p, method="BH")
  res <- cbind(id=rownames(et$table), et$table, adj.p, threshold=-p)
  rownames(res) <- NULL
  if (reorder){
    res <- res[order(res$id), ]
    message("Reordering genes!!!")
  }
  
  return(list(result=res, object=d, testobject=et))
}



run.deseq <- function(counts, labels,method="pooled",norm=FALSE, reorder=TRUE) {
  require(DESeq, quietly=TRUE)
  labels <- as.factor(labels)
  cds <- newCountDataSet(counts, labels)
  if (norm){
    cds <- rep(1,length(labels))
  } else {
    cds <- estimateSizeFactors(cds)
  }
  cds <- tryCatch(estimateDispersions(cds,method=method),
                  error=function(e) estimateDispersions(cds, fitType="local"))
  res <- nbinomTest(cds, levels(labels)[1], levels(labels)[2])
  colnames(res)[colnames(res) == "padj"] <- "adj.p"
  colnames(res)[colnames(res) == "log2FoldChange"] <- "logFC"
  res <- cbind(res, threshold=-res$pval)
  if (reorder){
    res <- res[order(res$id), ]
    message("Reordering genes!!!")
  }
  return(list(result=res, object=cds))
}


run.bayseq <- function(counts, labels, ncores=1, normmethod="quantile", reorder=TRUE) {
  require(baySeq, quietly=TRUE)
  require(parallel, quietly=TRUE)
  cl <- NULL
  if (ncores > 1) {
    cl <- makeForkCluster(nnodes=ncores)
    print(cl)
  }
  groups=list(nde=rep(1, ncol(counts)), de=as.numeric(as.factor(labels)))
  cd <- new("countData", data = counts, replicates = labels, groups = groups)
  
  if (normmethod=="none"){
    cd@libsizes <- rep(1,ncol(counts))
  } else {
    cd@libsizes <- getLibsizes(cd, estimationType=normmethod)
  }
  cd <- getPriors.NB(cd, samplesize = 1e5, estimation = "QL", cl = cl)
  cd <- getLikelihoods.NB(cd, pET = 'BIC', cl = cl)
  t <- topCounts(cd, number=nrow(counts), group = "de")
  colnames(t)[colnames(t) == "FDR"] <- "adj.p"
  res <- cbind(id=rownames(t), t, threshold=t$Likelihood)
  rownames(res) <- NULL
  if (ncores > 1)
    stopCluster(cl)
  if (reorder){
    res <- res[order(res$id), ]
    message("Reordering genes!!!")
  }	
  return(list(result=res, object=cd))
}



run.noiseq <- function(counts, labels, normmethod="tmm",reorder=TRUE) {
  if (!file.exists("noiseq.r"))
    system("wget http://bioinfo.cipf.es/noiseq/lib/exe/fetch.php?media=noiseq.r -O noiseq.r")
  source("noiseq.r")
  d1 <- as.matrix(counts[, labels == levels(labels)[1]])
  d2 <- as.matrix(counts[, labels == levels(labels)[2]])
  res <- noiseq(d1, d2, repl = "bio", k = 0.5, norm = normmethod, long = 1000, q = 0.8, nss = 0, lc = 1)
  p <- sort(res$probab, decreasing=TRUE,na.last=TRUE)
  fdr <- cumsum(1 - p)/1:length(p)
  res <- data.frame(id=names(p), adj.p = fdr, threshold=p)
  if (reorder){
    res <- res[order(res$id), ]
    message("Reordering genes!!!")
  }
  
  return(list(result=res))
}



run.degseq <- function(counts, labels, outputprefix="degseqtmp",reorder=TRUE) {
  suppressPackageStartupMessages(suppressWarnings(library(tcltk)))
  require(DEGseq, quietly=TRUE)
  dir <- tempfile(pattern=paste(outputprefix))
  m1 <- cbind(rownames(counts), counts[, labels == levels(labels)[1]])
  m2 <- cbind(rownames(counts), counts[, labels == levels(labels)[2]])
  DEGexp(method="MARS", outputDir=dir,
         geneExpMatrix1=m1, geneCol1=1, expCol1=2:ncol(m1), groupLabel1=levels(labels)[1],
         geneExpMatrix2=m2, geneCol2=1, expCol2=2:ncol(m2), groupLabel2=levels(labels)[2])
  res <- read.table(paste(dir, "output_score.txt", sep="/"), header=TRUE)
  unlink(dir, recursive=TRUE)
  colnames(res)[1] <- "id"
  colnames(res)[4] <- "logFC.unnorm"
  colnames(res)[5] <- "logFC"
  colnames(res)[8] <- "adj.p"
  colnames(res)[9] <- "adj.p2"
  colnames(res)[10] <- "is.de"
  res <- cbind(res, threshold=-log(res$p.value))
  if (reorder){
    res <- res[order(res$id), ]
    message("Reordering genes!!!")
  }
  
  return(list(result=res))
}


run.samseq <- function(counts, labels,reorder=TRUE) {
  library(samr, quietly=TRUE)
  x <- SAMseq(counts, labels, geneid=rownames(counts), resp.type="Two class unpaired", fdr.output=1.0)
  ss <- rbind(x$siggenes.table$genes.up, x$siggenes.table$genes.lo)
  res <- data.frame(id=ss[, 2], foldchange=as.numeric(ss[, 4]), q=as.numeric(ss[, 5]))
  res <- cbind(res, threshold=-res$q)
  missed.id = setdiff(rownames(counts), ss[, "Gene Name"])	
  res <- rbind(res, data.frame(id=missed.id, foldchange=0, q=max(res$q), threshold = 1-max(res$q)))
  res$id <- as.character(res$id)
  if (reorder){
    res <- res[order(res$id), ]
    message("Reordering genes!!!")
  }
  
  colnames(res)[colnames(res) == "q"] <- "adj.p"
  res$adj.p <- res$adj.p / 100
  return(list(result=res, object=x))
}



run.mixnb <- function(counts, labels=NULL, ncores=1,reorder=TRUE, ...) {
  library(mixnb)	
  object <- mixNB.parallel(counts, labels=labels, ...)
  if (is.null(labels) | missing("labels")){
    cat("Returning INI.\n")
    res <- data.frame(id=rownames(object$X), threshold=object$INI)
  } else{
    cat("Returning pval.\n")
    res <- data.frame(id=rownames(object$X), threshold=-object$pval, adj.p=object$pval)
  }
  if (reorder){
    res <- res[order(res$id), ]
    message("Reordering genes!!!")
  }
  res <- res[order(res$id), ]
  return (list(result=res, object=object))
}




run.dss <- function(counts, labels, normmethod="quantile",reorder=TRUE) {
  require(DSS, quietly=TRUE)
  idx <- which(rowSums(counts)>=1)
  X <- counts[idx, ]
  
  colnames(X) <- NULL
  d <- newSeqCountSet(X, as.numeric(labels)-1)
  if (normmethod != "none")
    d <- estNormFactors(d, method=normmethod)
  d <- estDispersion(d)
  object <- waldTest(d, 0, 1)
  
  id <- rownames(counts)
  
  threshold <- rep(0,nrow(counts))
  threshold[idx] <- abs(object$stats)[order(object$geneIndex)]
  
  res <- data.frame(id=id, threshold=threshold,stringsAsFactors=FALSE)
  if (reorder)
    res <- res[order(res$id), ]
  return(list(result=res,object=object))
}


run.poissonseq <- function(counts,labels,normmethod="default",
                           type="twoclass",reorder=TRUE){
  library(PoissonSeq)
  
  if (is.null(rownames(counts))){
    stop("Rownames of object \"counts\" must be set and unique!")
  }
  input <- list(n=counts,y=as.numeric(labels),pair=FALSE,type="twoclass",
                gname=rownames(counts))
  
  tmp <- tempfile(pattern="poissonseq-pow")
  para <- list(trans=FALSE,npermu=100,seed=10,ct.sum=5,ct.mean=0.5,div=10,pow.file=tmp)
  
  if (normmethod=="default")
    para$trans=TRUE
  
  obj <- PS.Main(input,para)
  
  pval <- rep(1,nrow(counts))
  names(pval) <- rownames(counts)
  pval[as.character(obj$gname)] <- as.numeric(obj$pval)
  
  res <- data.frame(id=rownames(counts), threshold=-log(pval))
  if (reorder)
    res <- res[order(res$id), ]
  
  return(list(result=res,object=obj))
}


run.deseq.mc <- function(counts, labels) {
  library(DESeq, quietly=TRUE, verbose=FALSE)
  design <- data.frame(row.names = colnames(counts), condition = labels)
  cds <- newCountDataSet(counts, design)
  cds <- estimateSizeFactors(cds) 
  cds <- estimateDispersions(cds, method="pooled-CR")
  fit1 <- fitNbinomGLMs(cds, count ~ condition)
  fit0 <- fitNbinomGLMs(cds, count ~ 1)
  pval <- nbinomGLMTest(fit1, fit0)
  
  res <- data.frame(id=rownames(fit1), threshold=-pval)
  return(list(result=res, object=cds, testobject=fit1))
}


run.edger.mc <- function(counts, labels, normmethod="none") {
  library(edgeR, quietly=TRUE)
  design <- model.matrix(~labels)
  d <- DGEList(counts=counts)
  if (normmethod != "none")
    d <- calcNormFactors(d, method=normmethod)
  d <- estimateGLMCommonDisp(d, design)
  d <- estimateGLMTrendedDisp(d,design)
  d <- estimateGLMTagwiseDisp(d, design)
  fit <- glmFit(d,design)
  lrt <- glmLRT(fit,coef=2:3) 
  
  p <- lrt$table$PValue
  res <- data.frame(id=rownames(lrt$table), threshold=as.numeric(-p))
  rownames(res) <- NULL
  return(list(result=res, object=d, testobject=lrt))
}


run.bayseq.mc <- function(counts, labels, normmethod="quantiles") {
  library(baySeq, quietly=TRUE)
  library(stringr)
  cl <- NULL
  labels <- as.numeric(as.factor(labels))
  groups=list(nde=rep(1, ncol(counts)), deAll=labels)
  for (i in unique(labels)) {
    tmp <- rep(1, length(labels))
    tmp[labels == i] <- 2
    n <- str_c("DE", i)
    groups[[n]] <- tmp
  }
  cd <- new("countData", data=counts, replicates=labels, groups=groups)
  if (normmethod=="none"){
    cd@libsizes <- rep(1,ncol(counts))
  } else {
    cd@libsizes <- getLibsizes(cd, estimationType=normmethod)
  }
  cd <- getPriors.NB(cd, samplesize = 1000, estimation = "QL", cl = cl)
  cd <- getLikelihoods.NB(cd, pET = 'BIC', cl = cl)
  t <- topCounts(cd, number=nrow(counts), group = "nde")
  res <- data.frame(id=rownames(t), threshold=as.numeric(-t$Likelihood))
  rownames(res) <- NULL
  return(list(result=res, object=cd, testobject=t))
}


run.samseq.mc <- function(counts, labels) {
  library(samr, quietly=TRUE)
  x <- SAMseq(counts, labels, geneid=1:nrow(counts), genenames=rownames(counts), resp.type="Multiclass", fdr.output=1.0)
  ss <- rbind(x$siggenes.table$genes.up, x$siggenes.table$genes.lo)
  res <- data.frame(id=ss[, 1], threshold=-as.numeric(ss[, 7])/100.0)
  missed.id = setdiff(rownames(counts), ss[, "Gene ID"])
  if (length(missed.id) > 0)
    res <- rbind(res, data.frame(id=missed.id, threshold = min(res$threshold)))
  return(list(result=res, object=x))
}


########## Dexus result to heatmap  #########################################
dexus2heatmap <- function(res, n=10, start=1, color="grey",
                          cexSamples=0.5,cexGenes=1,getGeneSymbols=FALSE,
                          newColNames=NULL,version=61,type="crosses",cexCrosses=2){
  library(RColorBrewer)
  
  idx <- order(res$INI,decreasing=TRUE)[(n+start-1):start]
  X <- log(res$X.norm[idx, ,drop=FALSE]+0.01)
  if (!is.null(newColNames)){
    colnames(X) <- newColNames
  }
  if (getGeneSymbols){
    rownames(X) <- getGeneSymbols(rownames(X),version=version)
  }
  
  
  m <- ncol(X)
  n <- nrow(X)
  resp <- t(res$resp[idx, ,drop=FALSE])
  
  X <- t(X)
  
  par(oma=c(1,6,0,0))
  
  rgb.palette <- colorRampPalette(c("white", "darkblue"),space = "rgb")
  image(X,xaxt="n",yaxt="n",col = rgb.palette(100))
  
  if (n==1){
    ia <- 1
    yy <- c(-1,1)
  } else {
    ia <- 1/(n-1)
    yy <- seq(-ia/2,1+ia/2,ia)
  }
  
  ib <- 1/(m-1)
  
  xx <- seq(-ib/2,1+ib/2,ib)
  
  
  for (j in 1:n){
    for (i in 1:m){
      
      if (type=="innerboxes"){
        if (resp[i,j]!=1){
          rect(xleft=xx[i]+ib/15,xright=xx[i+1]-ib/15,ybottom=yy[j]+ia/15,ytop=yy[j+1]-ia/15,
               lwd=5,border="red")
        }
      } else if (type=="boxes"){
        if (resp[i,j]!=1){
          rect(xleft=xx[i],xright=xx[i+1],ybottom=yy[j],ytop=yy[j+1],
               lwd=2,border="red")
        }
      } else {
        if (resp[i,j]!=1){
          points(x=(xx[i]+xx[i+1])/2,y=(yy[j]+yy[j+1])/2,pch=4,col="red",lwd=3,cex=cexCrosses)
        }
      }
    }	
  }
  
  axis(1,at=seq(0,1,ib),labels=rownames(X),las=2,cex.axis=cexSamples)
  if (n==1){axis(2,at=0,labels=colnames(X),las=2,cex.axis=cexGenes)  } else {axis(2,at=seq(0,1,ia),labels=colnames(X),las=2,cex.axis=cexGenes)
  } }

