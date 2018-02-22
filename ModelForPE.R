require(xlsx)

byapply <- function(x, by, fun, ...)
{
  # Create index list
  if (length(by) == 1)
  {
    nc <- ncol(x)
    split.index <- rep(1:ceiling(nc / by), each = by, length.out = nc)
  } else # 'by' is a vector of groups
  {
    nc <- length(by)
    split.index <- by
  }
  index.list <- split(seq(from = 1, to = nc), split.index)
  
  # Pass index list to fun using sapply() and return object
  sapply(index.list, function(i)
  {
    do.call(fun, list(x[, i], ...))
  })
}

df <- read.xlsx("Protein Table Trimester Category ALL the proteins.xlsx", sheetIndex = 1)
colnames(df)
rownames(df)<-df[,1] ;df<-df[,-c(1,2)]

ind <- apply(df, 1, var) == 0
df <- df[!ind,]
df<-(as.matrix(df))
y <- byapply(df, 3, rowMeans)

clinical <- read.xlsx("PREECLAMPSIA DATA FOR RAHUL.xlsx",sheetIndex = 1)
head(clinical)
colnames(clinical) <- gsub("\\.", "_", colnames(clinical))

table(clinical$X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia)

clinical$BMI_in_the_beginning_of_pregnancy<-as.numeric(as.character(
                        clinical$BMI_in_the_beginning_of_pregnancy))

clinical$BMI_at_the_end_of_pregnancy<-as.numeric(
              as.character(clinical$BMI_at_the_end_of_pregnancy))
w<-clinical$X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia

w[w==2]=1

clinical$X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia<-w

######################## Feature Selection###########################
library(relaimpo)
lmMod <- lm(X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia ~ 
              AGEATENTRY + Height + Weight_at_the_beginning_of_pregnancy+ 
              BMI_at_the_end_of_pregnancy +
              BMI_in_the_beginning_of_pregnancy , data = clinical)
relImportance <- calc.relimp(lmMod, type = "lmg", rela = TRUE)

lmMod <- lm(X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia ~ 
              AGEATENTRY + Height + Weight_at_the_beginning_of_pregnancy+ 
              BMI_in_the_beginning_of_pregnancy , 
            data = clinical)
selectedMod <- step(lmMod)
summary(lmMod)
#############################################################################

input_ones <- clinical[which(clinical$X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia
                              == 1), ]  # all 1's
input_zeros <- clinical[which(clinical$X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia
                               == 0), ]  # all 0's

set.seed(100)  # for repeatability of samples
input_ones_training_rows <- sample(1:nrow(input_ones), 1*nrow(input_ones))  # 1's for training
input_zeros_training_rows <- sample(1:nrow(input_zeros), 1*nrow(input_ones))  # 0's for training. Pick as many 0's as 1's
training_ones <- input_ones[input_ones_training_rows, ]  
training_zeros <- input_zeros[input_zeros_training_rows, ]
trainingData <-   rbind(training_ones, training_zeros)  # row bind the 1's and 0's 
dim(trainingData)
# Create Test Data
test_ones <- input_ones[sample(1:nrow(input_ones)), ]
test_zeros <- input_zeros[sample(1:nrow(input_zeros)), ]
testData <-   rbind(test_ones, test_zeros)  # row bind the 1's and 0's
view(testData)
library(PerformanceAnalytics)

chart.Correlation(clinical[,c(2,28,29,15,16)], 
                  method="spearman",
                  histogram=TRUE,
                  pch=16)

ind1= which(is.na(trainingData$BMI_at_the_end_of_pregnancy))
trainingData.omit <- trainingData[-ind1,]
ind2= which(is.na(testData$BMI_at_the_end_of_pregnancy))
testData.omit <- testData[-ind2,]

logitMod <- glm(X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia ~ 
                  BMI_in_the_beginning_of_pregnancy,
                  data=trainingData, family=binomial(link="logit"))

summary(logitMod)



predicted <- predict.glm(logitMod, testData, type="response")  # predicted scores
##predicted1 <- plogis(predict(logitMod, testData))
library(InformationValue)
optCutOff <- optimalCutoff(testData$X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia
                           , predicted)[1]
mis<-misClassError(testData$X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia
              , round(predicted), threshold = 0.01)
plotROC(testData$X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia
        , round(predicted))

con<-Concordance(testData$X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia
            , round(predicted))

sensitivity(testData$X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia, 
            round(predicted), threshold = 0.1)
#> 0.25
specificity(testData$X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia,
            round(predicted), threshold = 0.1)
#> 0.85

confusionMatrix(testData$X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia,
                round(predicted), threshold = optCutOff)

################################

#
ind_df<-grep("_3rd_",colnames(df))
y<-t(df[,ind_df])
y1<-cbind(y,clinical$X0_no_pre_eclampsia__1_Preeclampsia__not_severe__2__severe_preeclampsia)
y2<-cbind(y,clinical$BMI_in_the_beginning_of_pregnancy)
colnames(y1)[362]<-"PE"
colnames(y2)[362]<-"BMI"

##T-Test to correlated BMI with biomarker expression
g3<-y2[,362]
pv3 <- apply(y2[,-362], 2, function(x){ t.test(x ,g3 )$p.value})
adjp3<-p.adjust(pv3,method="hochberg")
mat3<-cbind(colnames(y2[,which(pv3<0.05)]),pv3[which(pv3<0.05)],adjp3[which(pv3<0.05)])
ord<-order(pv3[which(pv3<0.05)])
mat3<-mat3[ord,]
colnames(mat3)<-c("ID","pvalue","adj. pvalue")
write.table(mat3,paste("BMI_expression_ttest_","pval_3rd.csv",sep=""),sep=":")

##Wilcoxon.test to correlated PE with biomarker expression
wil<-wilcox.test(y1[,-362]~PE,data=y1,conf.int = TRUE)
g<-factor(as.numeric(y1[,362]))
pv <- apply(y1[,-362], 2, function(x){ wilcox.test(x ~g )$p.value})
adjp<-p.adjust(pv,method="hochberg")
mat<-cbind(colnames(y1[,which(pv<0.05)]),pv[which(pv<0.05)],adjp[which(pv<0.05)])
ord<-order(pv[which(pv<0.05)])
mat<-mat[ord,]
colnames(mat)<-c("ID","pvalue","adj. pvalue")

write.table(mat,paste("PE_expression_wilcox_","pval_3rd.csv",sep=""),sep=":")

require(coin)
pv.d <- apply(as.data.frame(y1[,-362]), 2, function(x){ pvalue(wilcox_test(x ~g,
                                                                        data=as.data.frame(y1[,-362]) ))})
adjp.d<-p.adjust(pv.d,method="hochberg")
mat.d<-cbind(colnames(y1[,which(pv.d<0.05)]),pv.d[which(pv.d<0.05)],adjp.d[which(pv.d<0.05)])
ord.d<-order(pv.d[which(pv.d<0.05)])
mat.d<-mat.d[ord.d,]
colnames(mat.d)<-c("ID","pvalue","adj. pvalue")  


##kruskal.test to correlated GROUP with biomarker expression
g1<-factor(rownames(y1))
pv1 <- apply(y1[,-362], 2, function(x){ kruskal.test(x ~g1 )$p.value})
adjp1<-p.adjust(pv1,method="hochberg")
mat1<-cbind(colnames(y1[,which(pv1<0.05)]),pv1[which(pv1<0.05)],adjp1[which(pv1<0.05)])
ord<-order(pv1[which(pv1<0.05)])
mat1<-mat1[ord,]
colnames(mat1)<-c("ID","pvalue","adj. pvalue")

write.table(mat1,paste("GROUP_expression_kruskal_","pval_3rd.csv",sep=""),sep=",")





