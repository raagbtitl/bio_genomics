## Code to run the network analysis for the PGC pathways paper
## Written by Neelroop Parikshak, neelroop@gmail.com
## R version 3.1.1 (2014-07-10) -- "Sock it to Me"
## Platform: x86_64-apple-darwin10.8.0 (64-bit)

## This code will allow one to reproduce the co-expression networks in O'Dushlaine et al., 2015

##### Goes through the following steps:
## 1) Get the dataset from GEO
## 2) Obtain updated annotations
## 3) Get necessary information for this study
## 4) Perform network analysis
## 5) Plot networks

## Obtain or load the necessary libraries

## Uncomment these lines if you do not have these packages
##source("http://bioconductor.org/biocLite.R")
##biocLite(c("WGCNA", "GEOquery", "igraph", "biomaRt"))

library(GEOquery) ## Automatically get GEO data
library(biomaRt) ## Annotations from biomaRt (host is set to use Gencode v10)
library(WGCNA) ## Weighted gene co-expresison network analysis package (v1.41.1)
library(igraph) ## Plot networks
options(stringsAsFactors=FALSE)

######## #1) Get data from Kang et al. on GEO
curr.study <- "KangEtAl"
gseid <- 25219
gsedat <- getGEO(paste("GSE",gseid,sep=""), GSEMatrix=TRUE)

datExpr <- exprs(gsedat[[1]]) ## Expression data
datMeta <- pData(phenoData(gsedat[[1]])) ## Phenotype data
geneProbeInfo <- pData(featureData(gsedat[[1]])) ## Probe information - let's use gene symbols

## Assign gene ID to each probe as given by the study's annotation - this could be improved in the future as it likely depends on Illumina's files from when the study was published
ENSGvec <- rep(NA,length=nrow(geneProbeInfo))
for (i in 1:length(ENSGvec)) {
  grout <- regexpr("gene:",as.character(geneProbeInfo[i,"mrna_assignment"]),fixed=TRUE)
  ENSGvec[i] <- substr(as.character(geneProbeInfo[i,"mrna_assignment"]),grout[1]+5,grout[1]+19)
}
keep <- substr(ENSGvec,1,4)=="ENSG"
datExpr <- datExpr[keep,]
geneDat <- ENSGvec[keep]
datMeta <- datMeta[,c(10:18)]

######## #2) Obtain updated annotations
## Use biomaRt to get Gencode v10 compatible probe information
listMarts(host="dec2011.archive.ensembl.org")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="dec2011.archive.ensembl.org")
identifier <- "ensembl_gene_id"
getinfo <- c("entrezgene","chromosome_name","start_position","end_position","ensembl_gene_id","hgnc_symbol","strand","gene_biotype")
geneDat <- getBM(attributes = getinfo,filters=identifier,values=ENSGvec[keep],mart=mart)
geneDat <- cbind(ENSGvec[keep],rownames(datExpr)[keep],geneDat[match(ENSGvec[keep],geneDat[,"ensembl_gene_id"]),])
rownames(datExpr) <- seq(1,nrow(datExpr),by=1)
rownames(geneDat) <- rownames(datExpr)

## Now collapse the rows - take the probe with the maximum mean to reflect the changes for any given gene
collapseDat <- collapseRows(datET = datExpr,rowGroup = geneDat[,"ENSGvec[keep]"],rowID = rownames(datExpr))
datExpr <- collapseDat$datETcollapsed
geneDat <- geneDat[as.numeric(collapseDat$group2row[,2]),]
match(colnames(datExpr),rownames(datMeta)) ## Samples match
subNames <- paste(substr(datMeta[,1],13,18),substr(datMeta[,2],9,12),substr(datMeta[,3],13,13),sep=".")
rownames(datMeta) <- colnames(datExpr) <- subNames
save(datExpr,datMeta,geneDat,file="./AllKangEtAllData_CollapsedRows_16874genes_1340samples_9metavars.Rdata")

######## #3) Get necessary information for this study
## Get the genes in the three pathways
immune <- c("ABCC8","ABL1","ACE","ACSL5","ACVR1","ACVR1B","ACVR1C","ACVR2B","ADAM10","ADAM9","ADCY1","ADCY2","ADCY3","ADCY4","ADCY6","ADCY8","ADCY9","ADIPOQ","AHSG","AKAP1","AKT1","AKT3","ALDOB","ANAPC1","ANAPC13","AP3S1","APBA1","APH1B","APOBEC1","APP","ARHGAP27","ARHGEF12","ATF2","ATP6V0A1","ATP6V0A4","ATP6V0D1","ATP6V0D2","ATP6V0E1","ATP6V1C1","ATP6V1C2","ATP6V1G1","ATP6V1H","BACE1","BAG1","BAG4","BAIAP2","BCAR1","BCL10","BCL2","BDKRB2","BDNF","BID","BLNK","BMP1","BMP2","BMP6","BMP7","BMP8A","BMP8B","BMPR1A","BMPR1B","BRAF","C14orf153","C3","CAB39L","CACNA1A","CACNA1B","CACNA1C","CACNA1D","CACNA1E","CACNA1G","CACNA1S","CACNB2","CALM2","CAMK2B","CAMK2D","CAMK2G","CAMK4","CAPN10","CASP10","CASP9","CCL5","CCND2","CCNE1","CD247","CD3D","CD3G","CD79B","CDAN1","CDC16","CDC25A","CDC42","CDH1","CDK6","CHRM2","CHRM3","CHRM5","CHRNA7","CHUK","CLDN10","CLDN16","CLDN20","CLEC4E","CPEB1","CR1","CREB1","CREB3L1","CREB3L2","CREB3L4","CREB5","CREBBP","CRK","CSK","CTNNB1","CTSS","DCLK1","DCP1A","DCP1B","DOK5","E2F3","EEA1","EGFR","EIF2AK2","EIF2AK3","EIF2AK4","EIF3E","EIF4E","ENG","ENPP1","EP300","ERBB2","ERBB3","ERBB4","ETS1","ETS2","ETV6","EXOC2","FCER1G","FCGR3A","FGF12","FGF14","FGF21","FGF22","FGF5","FGF6","FGFR1","FGFR2","FKBP1A","FKBP1B","FLRT2","FMN2","FOXO3","FRS2","FURIN","FZR1","GCK","GDNF","GHR","GNA11","GNA14","GNAI1","GNAI3","GNAL","GNAO1","GNAQ","GNAS","GNAZ","GPAM","GRB10","GRB14","GRB2","GSK3B","HBEGF","HDAC5","HDAC9","HGF","HKDC1","HLA-DMB","HLA-DOB","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DRA","HSP90AA1","HSP90AB1","HSP90B1","HSPA9","HSPD1","IFIT1","IFNAR1","IFNAR2","IFNGR2","IGF1","IGF1R","IKBKB","IKBKE","IL10","IL10RA","IL12B","IL18","IL1RL1","IL4R","IL6","INHBA","INHBC","INHBE","INPPL1","INSR","INSRR","IRS1","IRS2","ITGAX","ITGB2","ITPR1","ITPR2","ITPR3","JAK1","JAK2","JUN","KAT2B","KCNH8","KCNJ11","KDR","KIT","KL","KLC2","KLKB1","KRAS","KSR1","LAMP1","LDLR","LPIN1","LYN","MAD1L1","MAD2L2","MAP2K2","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K4","MAP3K5","MAP3K7","MAP3K9","MAP4K1","MAPK1","MAPK10","MAPK11","MAPK13","MAPK3","MAPK6","MAPK9","MDM2","MITF","MSTN","NCR1","NCSTN","NDUFA2","NFATC1","NFATC2","NFATC3","NFKB1","NFYA","NGF","NOD2","NOTCH1","NOTCH2","NOTCH3","NOTCH4","NRAS","NTF3","NTRK1","NTRK2","NTRK3","PAK1","PAK2","PARP1","PCSK1","PCSK2","PCSK5","PCSK6","PCSK7","PCSK9","PDE3A","PDE3B","PDGFB","PDGFC","PDGFD","PDGFRB","PGR","PHIP","PIAS2","PIGU","PIK3C3","PIK3CG","PIK3R1","PIP4K2B","PKMYT1","PKN1","PKN2","PLA2G12A","PLA2G2D","PLA2G2F","PLA2G3","PLA2G4E","PLCB1","PLCB2","PLCB4","PLCE1","PLCG1","PLCG2","PLD1","PLD2","PPARA","PPARG","PPAT","PPM1A","PPP2CB","PPP2R1A","PPP2R2A","PPP2R2B","PPP2R2C","PPP2R2D","PPP3CA","PPP3R1","PPP4C","PRKAA1","PRKAG1","PRKAG2","PRKCA","PRKCE","PRKCH","PRKCI","PRKCQ","PRKCZ","PSEN1","PTK2B","PTPN2","PTPRE","RAB5A","RAB5B","RAC1","RAF1","RALB","RAPGEF1","RGL1","RGS12","RHOA","RNASEL","RNMT","RPE65","RPS6","RPS6KA2","RPS6KA5","RPS6KB1","RRAS2","RXRA","SCARB1","SEL1L","SHC3","SHC4","SKI","SLC2A2","SLC2A4","SLC2A8","SMAD1","SMAD3","SMAD6","SMAD7","SMURF1","SNAI1","SORBS1","SORT1","SOS1","SOS2","SPIRE1","SPIRE2","SRD5A2","SRF","STK39","STXBP4","SYK","TBK1","TCF7","TCF7L1","TCF7L2","TGFA","TGFBR2","TIAM1","TIAM2","TIRAP","TLL1","TLL2","TLR2","TLR9","TNF","TP73","TRAIP","TSC1","TSC2","TYK2","TYRO3","UBE2B","UCP2","USF1","VAV1","VAV2","VAV3","VDR","YWHAE","YWHAG","YWHAH","YWHAZ","ZAK","ZFP106","ZFYVE16","ZFYVE9","ZNF274")

histone <- c("ACTL6A","ALDH1L2","ALKBH8","AMT","ARRB1","AS3MT","ASH1L","ASZ1","ATF7IP","BAZ2A","BHMT","BRCA2","BRD1","C13orf39","C17orf79","C2orf56","CARM1","CCDC101","COMT","CPA4","CREBBP","CSRP2BP","CTBP1","CXXC1","DFFB","DMGDH","DNMT3A","DNMT3L","DOHH","DOT1L","EHMT1","EHMT2","EIF5A2","ELP4","EP300","EP400","ETF1","FBXO11","FTCD","FTSJ2","FTSJ3","GART","GATAD2A","GNAS","GSTCD","GTF3C4","HAT1","HDAC4","HDAC9","HELLS","IRF4","JARID2","KAT2B","LCMT1","LCMT2","MAEL","MAP3K7","MAPK3","METT5D1","METTL4","METTL8","MGMT","MLL","MLL4","MLL5","MTFMT","MYOCD","MYST1","MYST2","MYST3","N6AMT2","NOS1","NSUN2","NSUN4","NSUN6","PAX5","PCMTD2","PEMT","PHF15","PHF17","PPARGC1A","PRDM2","PRDM5","PRMT2","PRMT5","PRMT7","PRMT8","RBBP5","RNMT","RPS6KA5","RUVBL2","SARDH","SETD1A","SETD6","SETD7","SETDB1","SETDB2","SF3B3","SHMT1","SMYD2","SMYD3","SNCA","SUPT3H","SUPT7L","SUV420H1","SUZ12","TAF5L","TAF6L","TARBP1","TCF3","TDRD5","TDRD9","TFB1M","TFB2M","TGS1","THUMPD2","THUMPD3","TRDMT1","TRRAP","VSIG2","WDR82","WHSC1","WHSC1L1","YEATS2","ZFP57")

synapse <- c("ABCA1","ACE","ACTG1","ACTN1","ACTN3","ACTR3","ACTR3B","ACVR1B","ADAM10","ADD1","ADORA1","ADORA2B","ADRA1A","ADRA1D","AGXT","AHSG","AIF1","AKT1","ALDOB","ALS2","ANKS1B","APBB2","APC","APP","AQP2","ASNS","ATG12","ATG16L1","ATG3","ATG9B","ATP2B4","ATP6V0E1","BCAR1","BCL11A","BCL2","BCL6","BDKRB2","BECN1","C12orf44","CACNA1A","CACNA1C","CADM1","CALD1","CAMK2D","CAPN5","CAPRIN2","CAPZA1","CAPZA2","CASR","CAV1","CCL24","CCNE1","CD38","CDA","CDC42EP2","CDH1","CDH10","CDH11","CDH12","CDH13","CDH15","CDH18","CDH2","CDH3","CDH4","CDH5","CDH7","CDH8","CDH9","CFTR","CGREF1","CHMP1A","CHPT1","CHRM3","CHRNA3","CLDN20","CNN3","CNTN2","COL17A1","CPEB1","CRB3","CREB1","CRIM1","CRLF3","CTNNA1","CTNNB1","CTNND1","CXADR","CYFIP1","DBH","DCC","DDN","DDR1","DERL2","DGKD","DLG1","DLG2","DLGAP1","DLGAP2","DLGAP3","DNM2","DNM3","DRD2","DSCAM","DST","DTNBP1","EAF2","ECE1","EDN1","EDNRA","EDNRB","EGLN2","EI24","EIF2AK4","ENO1","ENPP1","ENTPD5","EP300","EPB41L5","EPB49","EPHX2","ERBB2IP","ERBB4","ESR2","EXOC4","F11R","F2","F2R","FADS1","FBXO45","FGFR1","FGFR1OP","FGFR2","FGFRL1","FLNC","FOXN4","FYB","GABBR1","GIPC1","GJB2","GNG4","GRB2","GREM1","GRIA1","GRIA4","GRID2","GRIN2A","GRIN2B","GRIN3A","GRM1","GRM3","GRM5","GSN","HBEGF","HDAC9","HFE","HMOX1","HNF1B","HNF4A","HOMER2","HTR2A","HTRA1","ICAM1","IGF2BP1","IGFBP4","IGFBP7","IL17RB","IL7R","INADL","ING1","INHBA","ITCH","ITGA1","ITGA5","ITGA6","ITGB4","JUB","JUN","KANK1","KCND2","KDR","KEL","KIAA1109","KIAA1543","KIF26A","KIFC3","KNG1","LAMA4","LAMA5","LAMC1","LGI1","LIMA1","LIMS1","LIN7A","LIN7B","LMX1A","LPAR1","LPHN1","LRRK2","LTBP4","LYN","MACF1","MAD2L2","MAP1B","MAP1LC3A","MAP2K5","MAP3K1","MAPK8IP2","MAPT","MAX","MLL4","MLL5","MLLT4","MPP5","MPP7","MTDH","MTPN","MUL1","MYH10","MYL2","MYO5B","MYOCD","NCK2","NCOA2","NDEL1","NDRG4","NDUFA13","NEB","NEBL","NETO1","NGF","NLGN1","NLGN2","NOS1","NOS3","NOTCH2","NPC1","NPPA","NRCAM","NRG1","NRG3","NRP1","NTN1","NTRK3","NUDT1","NUMB","NUMBL","NUPR1","OMG","OSGIN1","P2RX4","P2RX6","PAFAH1B1","PAPPA2","PARD3","PARVA","PARVB","PCSK9","PDE5A","PDGFB","PDLIM5","PHF17","PIK3C3","PLCE1","PLEKHA7","PLXNA4","PMAIP1","PPARD","PPARG","PPP1R9B","PPP2R1A","PRKCI","PRKCQ","PRMT2","PSD3","PSEN1","PTGS2","PTK2B","PTPRK","PTPRM","PVRL1","PVRL3","PVRL4","RASA1","RERG","RICTOR","RRAGC","RSU1","RXRB","RYK","SCIN","SEMA3A","SERTAD2","SGMS1","SHANK1","SHANK2","SHANK3","SIPA1L1","SIRPA","SIRPB1","SIRPG","SKAP2","SLC1A2","SLC1A3","SLC8A1","SLIT1","SLIT2","SLIT3","SMAD3","SMAD7","SMARCA2","SMARCA4","SMO","SOD2","SORBS1","SPTA1","SPTBN1","SPTBN2","SQSTM1","SREBF1","SSH2","STRA8","STRN","TACR1","TANC1","TAOK2","TBCD","TBXAS1","TCHP","TGFBR3","TJP1","TLN1","TLN2","TNN","TP73","TRIM24","TRIP6","TRPV4","TSC1","TSC2","TSG101","USF1","USP47","VASP","VDR","VSIG2","WFDC1","WHSC1L1","WIPI1","WIPI2","WISP1","WISP3","WNT3","WNT3A","WNT5A","WNT6","WNT7A","WNT9B","YWHAZ","ZNF639") 

######## #4) Perform network analysis
## Load the expression data and store some metadata variables
load("./AllKangEtAllData_CollapsedRows_16874genes_1340samples_9metavars.Rdata")

sampstage <- as.numeric(substr(as.character(datMeta[,"characteristics_ch1.5"]),8,10)) ## Developmental stage

sampreg <- substr(as.character(datMeta[,"characteristics_ch1.1"]),9,20) ## Brain region

rmreg <- is.na(match(sampreg,names(table(sampreg))[table(sampreg)<20])) ## Remove regions that are too few in number
datExpr <- datExpr[,rmreg]
datMeta <- datMeta[rmreg,]

sampstage <- as.numeric(substr(as.character(datMeta[,"characteristics_ch1.5"]),8,10))
sampreg <- substr(as.character(datMeta[,"characteristics_ch1.1"]),9,20)

netVec <- vector(mode="list",length=1)
netVec[[1]]$list <- intersect(union(union(synapse,histone),immune),geneDat[,"hgnc_symbol"]) ## Intersection of the available genes with the union of the three broad pathways above

## Set an index to store the data
i=1
keep <- na.omit(match(netVec[[i]]$list,geneDat[,"hgnc_symbol"]))
netVec[[i]]$datExpr <- datExpr[keep,]
netVec[[i]]$geneDat <- geneDat[keep,]

## Now running WGCNA with signed bicor networks
powers <- seq(2,30,by=2)
sft <- pickSoftThreshold(t(netVec[[i]]$datExpr),
                         powerVector=powers,
                         corFnc="bicor",networkType="signed")

##   Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
##1      2   0.0377 -1.80          0.858 219.000  217.0000 255.00
##2      4   0.5930 -1.64          0.793  81.100   76.9000 128.00
##3      6   0.8260 -1.24          0.901  36.400   32.2000  78.20
##4      8   0.8890 -1.24          0.962  18.700   14.8000  54.70
##5     10   0.9410 -1.26          0.988  10.600    7.2300  40.30
##6     12   0.9430 -1.28          0.980   6.420    3.7500  30.80
##7     14   0.9240 -1.30          0.947   4.120    2.0300  24.10
##8     16   0.9270 -1.32          0.946   2.750    1.1100  19.20
##9     18   0.9530 -1.30          0.970   1.900    0.6290  15.50 <- chose 18 as the R^2 fit begins to peak here - however neighboringz values are very similar in resultant clusters
##10    20   0.9540 -1.33          0.963   1.350    0.3570  12.70
##11    22   0.9420 -1.33          0.944   0.984    0.2120  10.50
##12    24   0.9140 -1.33          0.907   0.730    0.1270   8.74
##13    26   0.9340 -1.33          0.926   0.551    0.0783   7.35
##14    28   0.9650 -1.30          0.959   0.422    0.0493   6.22
##15    30   0.8690 -1.34          0.832   0.327    0.0299   5.29

power <- 18 ## Power of 18 chosen here, though many other powers result in the same modules

## Run WGCNA using the semi-automated function
net <- blockwiseModules(t(netVec[[i]]$datExpr),power=power,
                        deepSplit=4,minModuleSize=10,minKMEtoStay=0,
                        mergeCutHeight=0.25,detectCutHeight=0.99995,
                        corType="bicor",networkType="signed",pamStage=FALSE,
                        verbose=3,saveTOMs=TRUE,maxBlockSize=10000) ## Main parameters are to cut the tree finely (deepSplit = 4), have a minimum module size of 10, and merge any modules with eigengene correlations > 0.75 to each other. This is a signed network using the biweight midcorrelation (a robust formulation of the correlation)

netVec[[i]]$TOMinfo <- net ## Store network analysis results

## Annotate the region names
regnames <- c("MFC","DFC","OFC","VFC", ## Frontal cortex
              "M1C","S1C","IPC", ## Motor, somatosensory, and association areas
              "A1C","STC","ITC", ## Auditory/temporal
              "V1C", ## Visual
              "HIP","AMY","STR","MD","CBC") ## Subcortical regions

## Plot dendrogram with module colors and trait correlations as well as individual module plots
pdf(file=paste("./GlobalNetworkPlots.pdf",sep=""),width=12,height=8)

## Make a numeric matrix to store metadata information
nummat <- matrix(0,nrow=10,ncol=length(sampstage))
nummat[1,] <- sampstage ## First row contains sample stage information
nummat[2,!is.na(match(sampreg,regnames[1:4]))] <- 1 ## Frontal cortex
nummat[3,!is.na(match(sampreg,regnames[5:7]))] <- 1 ## Motor/somatosensory/association cortex
nummat[4,!is.na(match(sampreg,regnames[8:10]))] <- 1 ## Auditory/temporal cortex
nummat[5,!is.na(match(sampreg,regnames[11:11]))] <- 1 ## Visual cortex

for (n in 6:10) {
  nummat[n,regnames[7+n-1]==sampreg] <- 1 ## Subcortical regions
}
rownames(nummat) <- c("Stage","Frontal","MSA","AT","Visual",regnames[12:16])
nummat <- t(nummat)

print(head(nummat)) ## The metadata matrix

## Now calculate gene correlations to each metadata variable and convert to colors for visualization purposes
netVec[[i]]$geneSignificance <- bicor(nummat,t(netVec[[i]]$datExpr))
colnames(netVec[[i]]$geneSignificance) <- rownames(netVec[[i]]$datExpr)
netVec[[i]]$geneSigColors <- t(numbers2colors(t(netVec[[i]]$geneSignificance),,signed=TRUE,lim=c(-1,1),naColor="black"))
colnames(netVec[[i]]$geneSigColors) <- rownames(netVec[[i]]$datExpr)
rownames(netVec[[i]]$geneSignificance) <- rownames(netVec[[i]]$geneSigColors) <- colnames(nummat)  

## Mark pathway membership information for each gene
immune.info <- !is.na(match(netVec[[i]]$geneDat[,"hgnc_symbol"],immune))
histone.info <- !is.na(match(netVec[[i]]$geneDat[,"hgnc_symbol"],histone))
synapse.info <- !is.na(match(netVec[[i]]$geneDat[,"hgnc_symbol"],synapse))

listLabels <- labels2colors(cbind(immune.info,histone.info,synapse.info))
colnames(listLabels) <- c("immune","histone","synapse")

## Plot a dendrogram with this information marked
plotDendroAndColors(dendro=netVec[[i]]$TOMinfo$dendrograms[[1]],
                    colors=t(rbind(netVec[[i]]$TOMinfo$colors,t(listLabels),netVec[[i]]$geneSigColors)),
                    cex.dendroLabels=1.2,addGuide=TRUE,
                    dendroLabels=FALSE,
                    groupLabels=c("Module Colors",colnames(listLabels),rownames(netVec[[i]]$geneSignificance)),
                    main=paste("Co-expression based topological overlap distance dendrogram",sep=" "))

## Plot eigengene dendrogram/heatmap - using bicor
MEs <- netVec[[i]]$TOMinfo$MEs
if (ncol(MEs) > 2) {
  plotEigengeneNetworks(MEs, "Eigengene (1st PC of each module) Clustering",
                        marHeatmap = c(3,4,2,2), marDendro = c(0,4,2,0),
                        plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))
}
colnames(MEs) <- substr(colnames(MEs),3,100)
rownames(MEs) <- colnames(netVec[[i]]$datExpr)


## Get sigend kME values - these delinate gene membership for each module in a continuous manner
tmpMEs <- MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
kMEdat <- signedKME(t(netVec[[i]]$datExpr), tmpMEs, corFnc="bicor")
netVec[[i]]$kMEtable <- kMEdat

######## #5) Plot networks
keepgenes <- NULL ## Store genes for summary plot in this variable

netTable <- cbind(netVec[[i]]$geneDat,netVec[[i]]$TOMinfo$colors,netVec[[i]]$kMEtable) ## Summarize network analysis results
netTable <- netTable[,-c(1:6,9,10)]
colnames(netTable)[3] <- "moduleColors"
write.table(cbind(netTable,immune.info,histone.info,synapse.info),file=paste("./kMEtable.txt",sep=""),sep=",") ## Store network analysis results

for (n in 1:(ncol(MEs)-1)) {
  ## Set plot layout
  layout(matrix(c(1,1,1,1,1,2,4,6,
                  1,1,1,1,1,2,4,6,
                  1,1,1,1,1,2,4,6,
                  1,1,1,1,1,3,5,7,
                  9,8,8,8,9,3,5,7,
                  9,8,8,8,9,3,5,7), ncol=6, nrow=8))
  keepcol <- netVec[[i]]$TOMinfo$colors == colnames(MEs)[n]
  
  ## Make an adjacency matrix and set some arbitrary cutoffs and parameters for plotting the network
  adjMat <- bicor(t(netVec[[i]]$datExpr[keepcol,]))
  maxsize <- min(25,nrow(adjMat))
  metsize <- min(10,nrow(adjMat))
  adjMat[adjMat<0.2] <- 0
  kME <- apply(adjMat,2,sum)
  cutoff <- sort(kME,decreasing=TRUE)[maxsize]
  cutoff2 <- sort(kME,decreasing=TRUE)[metsize]
  tophubs <- kME>=cutoff
  metgenes <- names(kME)[kME>=cutoff2]
  keepgenes <- c(keepgenes,metgenes)
  adjMat <- adjMat[tophubs,tophubs]
  numcors <- min(200,(maxsize^2-maxsize)/2)
  topcors <- sort(as.numeric(adjMat),decreasing=TRUE)[numcors]
  adjMat[adjMat<=topcors] <- 0
  
  ## Create a network object to store the data, plot these networks in a circle plot
  gA <- graph.adjacency(as.matrix(adjMat[1:5,1:5]),mode="undirected",weighted=TRUE,diag=FALSE) ## Top 5 in the center
  gB <- graph.adjacency(as.matrix(adjMat[6:maxsize,6:maxsize]),mode="undirected",weighted=TRUE,diag=FALSE) ## Additioanl genes on the periphery
  layoutCircle <- rbind(layout.circle(gA)/2,layout.circle(gB)) ## Construct layout
  
  g1 <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=TRUE,diag=FALSE) ## Graph information for plotting

  ## Color in the pie charts
  modColors <- t(listLabels)
  modColors[modColors=="turquoise"] <- "grey"
  modColors[1,modColors[1,]=="blue"] <- "darkblue"
  modColors[2,modColors[2,]=="blue"] <- "orange"
  modColors[3,modColors[3,]=="blue"] <- "black"

  valueList <- lapply(1:ncol(modColors), function(x) as.numeric(!is.na(modColors[,x])))
  colorList <- lapply(1:ncol(modColors), function(x) modColors[,x])

  ## Extract proper vertex names
  geneSymbols <- netVec[[i]]$geneDat[,"hgnc_symbol"]
  par(mar=c(1,1,1,1))

  ## Plot the network for this module, specifying plotting parameters
  plot.igraph(g1,vertex.label=geneSymbols[keepcol][tophubs],
              vertex.label.dist=0.8,
              vertex.shape="pie",
              vertex.pie=valueList,
              vertex.pie.color=colorList,
              vertex.size=12,
              vertex.label.color="black",
              vertex.color=colnames(MEs)[n],
              vertex.label.cex=1.2,
              vertex.frame.color="black",
              layout=layoutCircle,
              edge.color="green")
    
  par(mar=c(2,2,2,1))

  ## Plot MEs vs time by region
  yran <- range(MEs[,n])
  regvec <- colnames(nummat)[2:7]
  for (reg in regvec) {
    keepregs <- nummat[,reg]==1
    agenum <- nummat[keepregs,1]
    plotme <- MEs[keepregs,n]
    scatter.smooth(jitter(agenum),plotme,main=paste(reg,colnames(MEs)[n]),ylim=yran,xlim=c(0,16),
                   cex.main=0.7,degree=2,pch=19,cex=0.5,cex.lab=0.5,cex.axis=0.5)
  }
  par(mar=c(2,3,3,2))
  verboseBoxplot(x=MEs[,n],g=factor(sampreg,levels=regnames),main="ME difference by region",ylab="Module eigengene\n",las=2,xlab="")
  
}
dev.off()
  
## PLot top 10 genes from each module in a meta-network type layout
pdf("./NetworkSummaryGraph.pdf",width=12,height=8)
keepcol <- match(keepgenes,rownames(netVec[[i]]$datExpr))

adjMat <- adjMat2 <- bicor(t(netVec[[i]]$datExpr[keepcol,]))  
adjMat2[adjMat<0.3] <- 0

## Here we use a spring-based layout to highlight the modularity of the hub genes across modules
g1 <- graph.adjacency(as.matrix(adjMat2),mode="undirected",weighted=TRUE,diag=FALSE)  

set.seed(8675309) ## Set the random seed - this method is sensitive to the current state of the seed as it is randomly initialized
layoutFR <- layout.fruchterman.reingold(g1,dim=2) ## Spring based layout

## Color the pie nodes
modColors <- t(listLabels[keepcol,])
modColors[modColors=="turquoise"] <- "grey"
modColors[1,modColors[1,]=="blue"] <- "darkblue"
modColors[2,modColors[2,]=="blue"] <- "orange"
modColors[3,modColors[3,]=="blue"] <- "black"

valueList <- lapply(1:ncol(modColors), function(x) as.numeric(!is.na(modColors[,x])))
colorList <- lapply(1:ncol(modColors), function(x) modColors[,x])

colorVec <- netVec[[i]]$TOMinfo$colors[match(keepgenes,netVec[[i]]$geneDat[,"ensembl_gene_id"])]

## Extract proper vertex names
geneSymbols <- netVec[[i]]$geneDat[match(keepgenes,netVec[[i]]$geneDat[,"ensembl_gene_id"]),"hgnc_symbol"]
par(mar=c(1,1,1,1))

## Plot the network with specific parameters
plot.igraph(g1,vertex.label=geneSymbols,
            vertex.label.dist=0.3,
            vertex.shape="pie",
            vertex.pie=valueList,
            vertex.pie.color=colorList,
            vertex.size=5,
            vertex.label.color="black",
            vertex.frame.color=colorVec,
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR,
            edge.color="green")

dev.off()
