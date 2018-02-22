lambda <- 10^seq(10, -2, length = 100)
library(glmnet)
set.seed(489)
train<-data.matrix(trainingData)
ridge.mod <- glmnet(train[,c(2:14,28,15,16,18:27)], train[,17], alpha = 0, lambda = lambda)
#find the best lambda from our list via cross-validation
cv.out <- cv.glmnet(train[,c(2:14,28,15,16,18:27)], train[,17], alpha = 0)
bestlam <- cv.out$lambda.min
predict.glmnet(ridge.mod, s = 0, exact = T, type = 'coefficients')[1:6,]
#make predictions
ridge.pred <- predict(ridge.mod, s = bestlam, newx = data.matrix(testData[,c(2:14,28,15,16,18:27)]))

#check MSE
ytest<-data.matrix(testData[,17])
mean((round(ridge.pred)-round(ytest))^2)
aa<-predict(ridge.mod, type = "coefficients", s = bestlam)[1:15,]

lasso.mod <- glmnet(train[,c(2,28,15,16)], train[,17], alpha = 1, lambda = lambda)
lasso.pred <- predict(lasso.mod, s = bestlam, newx = data.matrix(testData[,c(2,28,15,16)]))
mean((round(lasso.pred)-ytest)^2)
lasso.coef  <- predict(lasso.mod, type = 'coefficients', s = bestlam)[1:5,]

kurtosis.test <- function (x) {
  m4 <- sum((x-mean(x))^4)/length(x)
  s4 <- var(x)^2
  kurt <- (m4/s4) - 3
  sek <- sqrt(24/length(x))
  totest <- kurt/sek
  pvalue <- pt(totest,(length(x)-1))
  pvalue 
}
kurtosis.test(df)
  
  skew.test <- function (x) {
    m3 <- sum((x-mean(x))^3)/length(x)
    s3 <- sqrt(var(x))^3
    skew <- m3/s3
    ses <- sqrt(6/length(x))
    totest <- skew/ses
    pt(totest,(length(x)-1))
    pval <- pt(totest,(length(x)-1))
    pval
  }
skew.test(df)