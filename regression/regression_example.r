# Copyright (C) 2019 Peter Sarvari, University of Southern California
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

library(glmnet)
df <- read.delim("13_methylation", header = TRUE, sep = "\t")
rownames(df) <- df$X
df$X <- NULL
df_s <- read.csv(file="filtered_successdatratio", header=TRUE, sep="\t")
df_s <- df_s[!(is.na(df_s$Female.Age) | df_s$Female.Age==""), ]
df_s <- df_s[!(is.na(df_s$Ratio) | df_s$Ratio==""), ]
rownames(df_s) <- df_s$Laboratory.ID
df_s$Laboratory.ID <- NULL

data <- merge(df, df_s, by=0)
rownames(data) <- data$Row.names
data$Row.names <- NULL

get.rsq <- function(y, y.hat) {
  sst <- sum((y - mean(y))^2)
  sse <- sum((y.hat - y)^2)
  return(1 - sse/sst)
}

cross.val.fold <- 10
n.individuals <- nrow(data)
print("Number of samples")
print(n.individuals)
parts <- sample(rep(1:cross.val.fold, length.out=n.individuals))

rsq_list <- rep(0,cross.val.fold)
varexp_train <- rep(0,cross.val.fold)
varexp_list <- rep(0,cross.val.fold)
mse_list <- rep(0,cross.val.fold)

for (i in 1:cross.val.fold) {

data$Numsuccess <- NULL
data$Numtrials <- NULL
training.individuals <- rownames(data)[which(parts != i)]
testing.individuals <- rownames(data)[which(parts == i)]
X_ <- data[training.individuals,]
Y <- as.matrix(X_$Ratio)
X_$Ratio <- NULL
X <- as.matrix(X_)

cv0 = cv.glmnet(X,Y,type.measure="mse",alpha=0,standardize=T)

if (i==1) {
    tmp_coeffs <- coef(cv0, s=cv0$lambda.min)
    df_coef <- data.frame(row.names = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

}
else {
    tmp_coeffs <- coef(cv0, s=cv0$lambda.min)
    df_coef2 <- data.frame(row.names = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
    df_coef <- merge(df_coef, df_coef2, by="row.names")
    rownames(df_coef) = df_coef$Row.names
    df_coef$Row.names <- NULL
    tmp <- paste("coefficient_", i)
    colnames(df_coef)[i] <- tmp
}

X.test <- data[testing.individuals,]
y_real <- as.matrix(X.test$Ratio)
X.test$Ratio <- NULL
x_pred <- as.matrix(X.test)

y_pred <- predict(cv0, newx = x_pred, s="lambda.min")
y.train.pred <- predict(cv0, newx=X, s="lambda.min")
result.lm = lm(y_real ~ y_pred)
print(cbind(y_real, y_pred))
rsq_list[i] <- summary(result.lm)$r.squared
mse_list[i] <- mean((y_real - y_pred)^2)
varexp_list[i] <- get.rsq(y_real, y_pred)
varexp_train[i] <- get.rsq(Y, y.train.pred)

}

df_coef
apply(df_coef, 1, mean)
print(paste("Mean R squared: ", mean(rsq_list)))
print(paste("Std dev of R squared: ", sd(rsq_list)))
print(paste("Mean MSE: ", mean(mse_list)))
print(paste("Std dev of MSE: ", sd(mse_list)))
print(paste("Mean of variance explained: ", mean(varexp_list)))
print(paste("Std dev of variance explained: ", sd(varexp_list)))
print(paste("Training set: Mean of variance explained: ", mean(varexp_train)))
print(paste("Training set: Std dev of variance explained: ", sd(varexp_train)))



