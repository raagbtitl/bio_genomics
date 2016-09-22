##########################
## Clustering Exercises ##
##########################

## Import a sample data set
## Download from GEO the Arabidopsis IAA treatment series "GSE1110" in TXT format. The direct link to the download is:
## ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/GSE1110/
## Uncompress the downloaded file.
## Import the data set into R
temp <- readLines("GSE1110_series_matrix.txt"); cat(temp[-grep("^!|^\"$", temp)], file="GSE1110clean.txt", sep="\n"); mydata <- read.delim("GSE1110clean.txt", header=T, sep="\t") # These import commands include a cleanup step to get rid of annotation lines and corrupted return signs.
rownames(mydata) <- mydata[,1]; mydata <- as.matrix(mydata[,-1]) # Assigns row indices and converts the data into a matrix object.

## Filtering
mydata <- mydata[apply(mydata>100, 1, sum)/length(mydata[1,])>0.5 & apply(log2(mydata), 1, IQR)>1.5,] # Retrieves all rows with high intensities (50% > 100) and high variability (IQR>1.5).

## Hierarchical clustering routine
y <- mydata
hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete"); hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete") # Generates row and column dendrograms.
mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] # Cuts the tree and creates color vector for clusters.
library(gplots); myheatcol <- redgreen(75) # Assign your favorite heatmap color scheme. Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]). Type demo.col(20) to see more color schemes.
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=myheatcol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc) # Creates heatmap for entire data set where the obtained clusters are indicated in the color bar.
x11(height=6, width=2); names(mycolhc) <- names(mycl); barplot(rep(10, max(mycl)), col=unique(mycolhc[hr$labels[hr$order]]), horiz=T, names=unique(mycl[hr$order])) # Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'.
clid <- c(1,2); ysub <- y[names(mycl[mycl%in%clid]),]; hrsub <- hclust(as.dist(1-cor(t(ysub), method="pearson")), method="complete") # Select sub-cluster number (here: clid=c(1,2)) and generate corresponding dendrogram.
x11(); heatmap.2(ysub, Rowv=as.dendrogram(hrsub), Colv=as.dendrogram(hc), col=myheatcol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc[mycl%in%clid]) # Create heatmap for chosen sub-cluster.
data.frame(Labels=rev(hrsub$labels[hrsub$order])) # Print out row labels in same order as shown in the heatmap.

## Hierarchical clustering
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R") # Import an alternative color scheme for the heatmap function.
mydatascale <- t(scale(t(mydata))) # Centers and scales data.
hr <- hclust(as.dist(1-cor(t(mydatascale), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete") # Clusters columns by Spearman correlation.
heatmap(mydata, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=my.colorFct(), scale="row") # Plot the data table as heatmap and the cluster results as dendrograms.
mycl <- cutree(hr, h=max(hr$height)/1.5) # Cut the tree at specific height and color the corresponding clusters in the heatmap color bar.
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
heatmap(mydata, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", RowSideColors=mycolhc) 

## Obtain significant clusters by pvclust bootstrap analysis
library(pvclust) # Loads the required pvclust package.
pv <- pvclust(scale(t(mydata)), method.dist="correlation", method.hclust="complete", nboot=10) # Perform the hierarchical cluster analysis. Due to time resrictions, we are using here only 10 bootstrap repetitions. Usually, one should use at least 1000 repetitions.
plot(pv, hang=-1); pvrect(pv, alpha=0.95) # Plots result as a dendrogram where the significant clusters are highlighted with red rectangles.
clsig <- unlist(pvpick(pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters) # Retrieve members of significant clusters.
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/dendroCol.R") # Import tree coloring function.
dend_colored <- dendrapply(as.dendrogram(pv$hclust), dendroCol, keys=clsig, xPar="edgePar", bgr="black", fgr="red", pch=20) # Create dendrogram object where the significant clusters are labeled in red.
heatmap(mydata, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", RowSideColors=mycolhc) # Plot the heatmap from above, but with the significant clusters in red and the cluster bins from the tree cutting step in the color bar.
x11(height=12); heatmap.2(mydata, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", trace="none", RowSideColors=mycolhc) # Plot heatmap with heatmap.2() function which scales better for many entries.
mydatasort <- mydata[pv$hclust$labels[pv$hclust$order], hc$labels[hc$order]] # Sort rows in data table by 'dend_colored' and its colums by 'hc'.
x11(height=16, width=12); par(mfrow=c(1,2)); plot(dend_colored, horiz=T, yaxt="n"); image(scale(t(mydatasort)), col=my.colorFct(), xaxt="n",yaxt="n") # Plot heatmap with bootstrap tree in larger format using instead of heatmap the image function.
pdf("pvclust.pdf", height=21, width=6); plot(dend_colored, horiz=T, yaxt="n"); dev.off(); pdf("heatmap.pdf", height=20, width=6); image(scale(t(mydatasort)), col=my.colorFct(), xaxt="n",yaxt="n"); dev.off() # Save graphical results to two PDF files: 'pvclust.pdf' and'heatmap.pdf'.

## Compare PAM (K-means) with hierarchical clustering
library(cluster) # Loads required library.
mydist <- t(scale(t(mydata))) # Center and scale data.
mydist <- as.dist(1-cor(t(mydist), method="pearson")) # Generates distance matrix using Pearson correlation as distance method.
pamy <- pam(mydist, max(mycl)) # Clusters distance matrix into as many clusters as obtained by tree cutting step (6).
mycolkm <- sample(rainbow(256)); mycolkm <- mycolkm[as.vector(pamy$clustering)]; heatmap(mydata, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", RowSideColors=mycolkm) # Compare PAM clustering results with hierarchical clustering by labeling it in heatmap color bar.
pdf("pam.pdf", height=20, width=20); heatmap(mydata, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", RowSideColors=mycolkm); dev.off() # Save graphical results to PDF file: 'pvclust.pdf'.

## Compare SOM with hierarchical clustering
library(som) # Loads required library.
y <- t(scale(t(mydata))) # Center and scale data.
y.som <- som(y, xdim = 2, ydim = 3, topol = "hexa", neigh = "gaussian") # Performs SOM clustering.
plot(y.som) # Plots results.
pdf("som.pdf"); plot(y.som); dev.off() # Save plot to PDF: 'som.pdf'.
somclid <- as.numeric(paste(y.som$visual[,1], y.som$visual[,2], sep=""))+1 # Returns SOM cluster assignment in order of input data.
mycolsom <- sample(rainbow(256)); mycolsom <- mycolsom[somclid]; heatmap(mydata, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", RowSideColors=mycolsom) # Compare SAM clustering results with hierarchical clustering by labeling it in heatmap color bar.
pdf("somhc.pdf", height=20, width=20); heatmap(mydata, Rowv=dend_colored, Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", RowSideColors=mycolsom); dev.off() # Save graphical results to PDF file: 'somhc.pdf'.

## Compare PCA with SOM
pca <- prcomp(mydata, scale=T) # Performs principal component analysis after scaling the data.
summary(pca) # Prints variance summary for all principal components.
library(scatterplot3d) # Loads 3D library.
scatterplot3d(pca$x[,1:3], pch=20, color=mycolsom) # Plots PCA result in 3D. The SOM clusters are highlighted in their color.
pdf("pcasom.pdf"); scatterplot3d(pca$x[,1:3], pch=20, color=mycolsom); dev.off() # Saves PCA plot in PDF format: 'pcasom.pdf'.

## Compare MDS with HC, SOM and K-means
loc <- cmdscale(mydist, k = 3) # Performs MDS analysis and returns results for three dimensions.
x11(height=8, width=8, pointsize=12); par(mfrow=c(2,2)) # Sets plotting parameters.
plot(loc[,1:2], pch=20, col=mycolsom, main="MDS vs SOM 2D") # Plots MDS-SOM comparison in 2D. The SOM clusters are highlighted in their color.
scatterplot3d(loc, pch=20, color=mycolsom, main="MDS vs SOM 3D") # Plots MDS-SOM comparison in 3D.
scatterplot3d(loc, pch=20, color=mycolhc, main="MDS vs HC 3D") # Plots MDS-HC comparison.
scatterplot3d(loc, pch=20, color=mycolkm, main="MDS vs KM 3D") # Plots MDS-KM comparison.

## Fuzzy clustering
library(cluster) # Loads cluster library.
fannyy <- fanny(mydist, k= max(mycl), memb.exp = 1.5); round(fannyy$membership, 2); fannyy$clustering # Performs fuzzy clustering with as many coefficients as clusters were obtained by tree cutting step in HC. The hard clustering results are provided in the 'clustering' slot.
fannyyMA <- round(fannyy$membership, 2) > 0.3; apply(fannyyMA, 1, which) # Returns multiple cluster memberships for coefficient above a certain value (here >0.3).



