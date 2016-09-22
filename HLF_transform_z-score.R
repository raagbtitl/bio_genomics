byapply <- function(x, by, fun, ...)
{  # Create index list
  if (length(by) == 1)
  {    nc <- ncol(x)
    split.index <- rep(1:ceiling(nc / by), each = by, length.out = nc)
  } else # 'by' is a vector of groups
  {    nc <- length(by); split.index <- by  }
  index.list <- split(seq(from = 1, to = nc), split.index)  
  # Pass index list to fun using sapply() and return object
  sapply(index.list, function(i)
  {    do.call(fun, list(x[, i], ...))  })}

a<-read.table("annotated_all_count_slimmed_nodup_upper_quartile.csv",header=T,sep="\t",row.names=1)
a<-as.matrix(a)

Cell_type<-c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3),rep(7,3),rep(8,3))
Cell_type<-factor(Cell_type)

design<-model.matrix(~1+Cell_type)
v<-voom(a,design,plot=TRUE)
hist(v$E)
v1<-byapply(v$E, 3, rowMeans)

write.table(v1,"mean_uq_1.csv",sep="\t")
v2<-scale(v1,center = TRUE, scale = TRUE) ##z-score
write.table(v2,"mean_uq_1_zscore.csv",sep="\t")



