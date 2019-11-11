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

# Rscript regression_v2.R ../input/AgingFull_Human_Sperm_hmr.txt
# ../input/Human_Sperm_hmr.txt ../input/AgingFullMeta.rds 
# ../input/agedat.tsv 1

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
stop('regression.R <methylation-data> <independent-data> 
     <age-data> <independent-age-data> <alpha>',
     call.=FALSE)
}
methylation.filename <- args[1]
newdata.filename <- args[2]
metadata.filename <- args[3]
newmetadata.filename <- args[4]
alpha.lasso <- args[5]

read.meta <- function(metadata.filename){
if (length(grep(".rds", metadata.filename))>0){
df.m <- readRDS(metadata.filename)
} else {
df.m <- read.csv(file=metadata.filename, header=TRUE, sep="\t")
}
if ("Laboratory.ID" %in% colnames(df.m)) {
rownames(df.m) <- df.m$Laboratory.ID
df.m$Laboratory.ID <- NULL
df.m$Age <- df.m$Male.Age
df.m$Male.Age <- NULL
}
return (df.m)
}

library(glmnet)
df <- read.table(methylation.filename, header=T)
id.df <- read.table(newdata.filename, header=T)
rownames(df) <- gsub("X", "", rownames(df))
rownames(id.df) <- gsub("X", "", rownames(id.df))
df.m <- read.meta(metadata.filename)
id.df.m <- read.meta(newmetadata.filename)

data <- merge(df, df.m, by=0)
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
if (alpha.lasso > 0) {
used.feat <- list()
}

#print("Check0")

for (i in 1:cross.val.fold) {

#print("Check01")

training.individuals <- rownames(data)[which(parts != i)]
testing.individuals <- rownames(data)[which(parts == i)]
#print("Check02")
X_ <- data[training.individuals,]
Y <- as.matrix(X_$Age)
X_$Age <- NULL
X <- as.matrix(X_)
#print("check03")

cv0 = cv.glmnet(X,Y,type.measure="mse",alpha=alpha.lasso,
                standardize=T)
#print("Check1")
if (alpha.lasso > 0) {
cs.k <- coef(cv0, s='lambda.min')
#print("Check2")
used.feat.k <- cs.k@Dimnames[[1]][cs.k@i + 1]
print(used.feat.k)
#print("Check3")
used.feat[[i]] <- used.feat.k
#print("Check4")
}

X.test <- data[testing.individuals,]
y_real <- as.matrix(X.test$Age)
X.test$Age <- NULL
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

print(used.feat)
used.feat.unlist <- unlist(used.feat)
used.unique <- unique(used.feat.unlist)
selected.features <- used.unique[-1] 
#excluding intercept
print(used.unique)
print(length(used.unique))

print(paste("Mean R squared: ", mean(rsq_list)))
print(paste("Std dev of R squared: ", sd(rsq_list)))
print(paste("Mean MSE: ", mean(mse_list)))
print(paste("Std dev of MSE: ", sd(mse_list)))
print(paste("Mean of variance explained: ", mean(varexp_list)))
print(paste("Std dev of variance explained: ", sd(varexp_list)))
print(paste("Training set: Mean of variance explained: ", 
            mean(varexp_train)))
print(paste("Training set: Std dev of variance explained: ", 
            sd(varexp_train)))
print(paste("Number of features: ", ncol(data)-1))


print("Training on second dataset using selected features")

id.data <- merge(id.df, id.df.m, by=0)
rownames(id.data) <- id.data$Row.names
id.data$Row.names <- NULL

featlen <- length(selected.features)
selected.features[featlen+1] <- "Age"
id.data <- id.data[, colnames(id.data) %in% selected.features]

n.individuals <- nrow(id.data)
print("Number of samples")
print(n.individuals)
parts <- sample(rep(1:cross.val.fold, length.out=n.individuals))

rsq_list <- rep(0,cross.val.fold)
varexp_train <- rep(0,cross.val.fold)
varexp_list <- rep(0,cross.val.fold)
mse_list <- rep(0,cross.val.fold)

for (i in 1:cross.val.fold) {

training.individuals <- rownames(id.data)[which(parts != i)]
testing.individuals <- rownames(id.data)[which(parts == i)]
X_ <- id.data[training.individuals,]
Y <- as.matrix(X_$Age)
X_$Age <- NULL
X <- as.matrix(X_)

cv0 = cv.glmnet(X,Y,type.measure="mse",alpha=1,
                standardize=T)

X.test <- id.data[testing.individuals,]
y_real <- as.matrix(X.test$Age)
X.test$Age <- NULL
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

print(paste("Found features:", ncol(id.data)))
print(paste("Total features:", length(selected.features)))

print(paste("Mean R squared: ", mean(rsq_list)))
print(paste("Std dev of R squared: ", sd(rsq_list)))
print(paste("Mean MSE: ", mean(mse_list)))
print(paste("Std dev of MSE: ", sd(mse_list)))
print(paste("Mean of variance explained: ", mean(varexp_list)))
print(paste("Std dev of variance explained: ", sd(varexp_list)))
print(paste("Training set: Mean of variance explained: ", 
            mean(varexp_train)))
print(paste("Training set: Std dev of variance explained: ", 
            sd(varexp_train)))
print(paste("Number of features: ", ncol(id.data)-1))

