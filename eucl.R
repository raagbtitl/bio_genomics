sq.Euc.dist <- function(x,y) {
  x <- as.matrix(x);  y <- as.matrix(y);  nr=nrow(x);  nc=nrow(y)
  x2 <- rowSums(x^2);  xsq = matrix(x2,nrow=nr,ncol=nc)
  y2 <- rowSums(y^2);  ysq = matrix(y2,nrow=nr,ncol=nc,byrow=TRUE)
  xy = x %*% t(y);  d = xsq + ysq - 2*xy
  if(identical(x,y)) diag(d) = 0;  d[which(d < 0)] = 0;  return(d)
}
data <-read.table("Filtered.csv",header=T,sep="\t");attach(data)
x<-sq.Euc.dist(Mo2,Mo3);write.table(x,"Mo2_Mo3_eucl.txt",sep="\t")
library(fields);x1<-rdist(Mo2,Mo3);write.table(x,"Mo2_Mo3_eucli.txt",sep="\t")

dist(rbind(data[,2], data[,3]))
