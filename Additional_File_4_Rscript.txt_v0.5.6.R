# script to reproduce graphics for the paper
# tested with R 2.15.2 and methylKit 0.5.6

# load the package
library(methylKit)


# data file is from reads aligned to hg18 assembly
# it contains per base % methylation values, coverage and strand information in addition to extra columns such as identifiers
download.file("http://methylkit.googlecode.com/files/bcERmyobj.rda",destfile="bcERmyobj.rda") # downloads the R data file, ?download.file for details or download manually
load("bcERmyobj.rda")

# FIGURE 2A
# plot methylation statistics on 2nd sample "T47D" in myobj which is a class of methylRawList
getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)


# FIGURE 2B
# plot coverage statistics on 2nd sample "T47D" in myobj which is a class of methylRawList
getCoverageStats(myobj[[2]],plot=TRUE,both.strands=FALSE)

# merge all samples to one table by using base-pair locations that are covered in all samples
# setting destrand=TRUE, will merge reads on both strans of a CpG dinucleotide. This provides better 
# coverage, but only advised when looking at CpG methylation
meth=methylKit::unite(myobj,destrand=FALSE)

# FIGURE 3
# plots pair-wise scatter plots from sample set
getCorrelation(meth,plot=TRUE)


# FIGURE 4A
# cluster all samples using correlation distance and plot hiarachical clustering
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

# FIGURE 4B
# principal component anlaysis of all samples.
# it plots loadings on PC1 and PC2 for each sample
# the samples that have similar loadings should be similar to each other
PCASamples(meth,scale=FALSE,center=FALSE)


#############################################################################
# EXTRA example: cluster using k-means clustering
#############################################################################
cl=kmeans(t(percMethylation(meth)),centers=2)
#-----------------------------------------------------------------------------



# calculate differential methylation p-values and q-values
myDiff=calculateDiffMeth(meth)

#############################################################################
# EXTRA example: calculate differential methylation with moderated t-test from limma package
#############################################################################
group<-factor(rep(c(0,1),each=2)) 
design<-model.matrix(~group)
design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,0,1,1,1))

# do the test in limma
library(limma)
p.meth=percMethylation(meth) # get percent methylation values
fit <- lmFit(p.meth, design = design)
fit2 <- ebayes(fit)

# make the data for methylKit object
df=cbind(meth[,1:5],pvalue=fit2$p.value[,2],qvalue=p.adjust(fit2$p.value[,2],method="BH"),meth.diff=rowMeans(p.meth[,1:2])-rowMeans(p.meth[,3:4])  )

# create a new methylDiff object
obj=new("methylDiff",df,sample.ids=meth@sample.ids,assembly=meth@assembly,context=meth@context,
        treatment=meth@treatment,destranded=meth@destranded,resolution=meth@resolution)
#------------------------------------------------------------------------------------


# FIGURE 5A
# plot percentages of differentially methylated bases over all chromosomes
diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=0.01,meth.cutoff=25,exclude=c("chrM","chrY","chrX") )


# get differentially methylated regions with 25% difference and qvalue<0.01
myDiff25p=get.methylDiff(myDiff,difference=25,qvalue=0.01 )

# FIGURE 5B
# get bedgraph file for differentially methylated bases
# ready to upload to UCSC
bedgraph(myDiff25p,file.name="diff.bedgraph",col.name="meth.diff")


# read-in the transcript locations to be used in annotation, a text file in BED format
# either full path to the text file or an URL to the BED file must be provided
gene.obj=read.transcript.features("http://methylkit.googlecode.com/files/refseq.hg18.bed.txt")


# FIGURE 6B
# annotate differentially methylated Cs with promoter/exon/intron using annotation data
gene.ann=annotate.WithGenicParts(myDiff25p,gene.obj)
plotTargetAnnotation(gene.ann)# plot pie chart for annotation



# FIGURE 6A
# plot nearest distance to TSS
tss.assoc=getAssociationWithTSS(gene.ann)
# plot the distance as histogram, plot only the ones that are at most 100kb away
hist(tss.assoc$dist.to.feature[abs(tss.assoc$dist.to.feature)<=100000],main="distance to nearest TSS",xlab="distance in bp",breaks=50,col="brown4")


# FIGURE 6C
cpgi.obj=read.feature.flank("http://methylkit.googlecode.com/files/refseq.hg18.bed.txt") # read the annotation file for CpG islands and shores, shores are determined by 2kb flanks
cpgi.ann=annotate.WithFeature.Flank(myDiff25p,cpg.obj[[1]],cpg.obj[[2]],feature.name="CpGi",flank.name="shores") # annotate differentially methylated bases with CpG islands and shores
plotTargetAnnotation(cpgi.ann) # plot pie chart for annotation

# FIGURE 6D
enhancer.obj=read.bed("http://methylkit.googlecode.com/files/all.enhancers.hg18.hmm.bed")  # read the annotation file for enhancers
enhancer.ann=annotate.WithFeature(myDiff25p,enhancer.obj,feature.name="enhancers") #  annotate differentially methylated bases with enhancers
plotTargetAnnotation(enhancer.ann) # plot pie chart for annotation

 
