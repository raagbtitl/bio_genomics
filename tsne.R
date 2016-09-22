library(tsne)

# initialize counter to 0
x <- 0
epc <- function(x) {
  x <<- x + 1
  filename <- paste("plot", x, "jpg", sep=".")
  cat("> Plotting TSNE to ", filename, " ")  
  # plot to file of 2400x1800 dimension
  jpeg(filename, width=2400, height=1800)  
  plot(x, t='n', main="T-SNE")
  text(x, labels=rownames(mydata))
  dev.off()
}
tsne_data <- tsne(d, k=5, epoch_callback=epc, max_iter=500, epoch=100)