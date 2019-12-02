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

id.data <- merge(id.df, id.df.m, by=0)
rownames(id.data) <- id.data$Row.names
id.data$Row.names <- NULL

col.names.data <- colnames(data)
col.names.id.data <- colnames(id.data)
col.names.same <- intersect(col.names.data, col.names.id.data)
id.data <- id.data[, colnames(id.data) %in% col.names.same]
data <- data[, colnames(data) %in% col.names.same]

print("Number of features in D1")
print(ncol(data))
print("Number of features in D2")
print(ncol(id.data))

get.rsq <- function(y, y.hat) {
  sst <- sum((y - mean(y))^2)
  sse <- sum((y.hat - y)^2)
  return(1 - sse/sst)
}

cross.val.fold <- 10
n.individuals <- nrow(data)
print("Number of samples")
print(n.individuals)
set.seed(1)
foldid <- sample(rep(1:cross.val.fold, length.out=n.individuals))

if (alpha.lasso > 0) {
used.feat <- list()
}

X_ <- data
Y <- as.matrix(X_$Age)

#Age centering
Y.mean <- mean(Y)
Y <- Y-Y.mean
print(Y.mean)

X_$Age <- NULL
X <- as.matrix(X_)
#Centering
X <- scale(X, center=TRUE, scale=FALSE)

cv0 = cv.glmnet(X,Y,type.measure="mse",alpha=alpha.lasso,
                nfolds=cross.val.fold, foldid = foldid)
if (alpha.lasso > 0) {
cs.k <- coef(cv0, s='lambda.min')
used.feat <- cs.k@Dimnames[[1]][cs.k@i + 1]
print(used.feat)
}

X.test <- id.data
y_real <- as.matrix(id.data$Age)

#Calculating new data mean; if we use this, then we use the unknown Y's
#But in a clinical setting they're obviously known
y_real.mean <- mean(y_real)

X.test$Age <- NULL
x_pred <- as.matrix(X.test)
#Centering
x_pred <- scale(x_pred, center=TRUE, scale=FALSE)

print("Predicting on second dataset using coefficients and best lambda")

y_pred <- predict(cv0, newx = x_pred, s="lambda.min")
y.train.pred <- predict(cv0, newx=X, s="lambda.min")

#Age add mean back
y_pred <- y_pred + y_real.mean
y.train.pred <- y.train.pred + Y.mean
Y <- Y + Y.mean

result.lm = lm(y_pred ~ y_real)
print(cbind(y_real, y_pred))

figure.filename <- sprintf('%s_agecenter.pdf', newdata.filename)

rsq_list <- summary(result.lm)$r.squared
mse_list <- mean((y_real - y_pred)^2)
varexp_list <- get.rsq(y_real, y_pred)
varexp_train <- get.rsq(Y, y.train.pred)

pdf(figure.filename, width=6.5, height=8, paper='US')
template <- "Independent data: \n %s \n (R2=%.3f)"
plot(y_real, y_pred, xlab="Actual age", ylab="Predicted age",
     main=sprintf(template, newdata.filename, varexp_list))
abline(result.lm)
abline(a=0, b=1, col="red", lty=2)
dev.off()

print(paste("Number of features selected by Lasso: ",
            length(used.feat)-1)) #-1 because of intercept


print(paste("Mean R squared: ", mean(rsq_list)))
print(paste("Mean MSE: ", mean(mse_list)))
print(paste("Mean of variance explained: ", mean(varexp_list)))
print(paste("Training set: Mean of variance explained: ", 
            mean(varexp_train)))
print(paste("Number of features: ", ncol(data)-1))

print("Training on second dataset using selected features")

selected.features <- used.feat
featlen <- length(selected.features)
selected.features[featlen+1] <- "Age"
id.data <- id.data[, colnames(id.data) %in% selected.features]

n.individuals <- nrow(id.data)
print("Number of samples")
print(n.individuals)
set.seed(1)
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

n.training <- nrow(X)
set.seed(1)
foldid <- sample(rep(1:cross.val.fold, length.out=n.training))

cv0 = cv.glmnet(X,Y,type.measure="mse",alpha=alpha.lasso,
                nfolds = cross.val.fold, foldid = foldid)

X.test <- id.data[testing.individuals,]
y_real <- as.matrix(X.test$Age)
X.test$Age <- NULL
x_pred <- as.matrix(X.test)

y_pred <- predict(cv0, newx = x_pred, s="lambda.min")
y.train.pred <- predict(cv0, newx=X, s="lambda.min")
result.lm = lm(y_real ~ y_pred)
#print(cbind(y_real, y_pred))
rsq_list[i] <- summary(result.lm)$r.squared
mse_list[i] <- mean((y_real - y_pred)^2)
varexp_list[i] <- get.rsq(y_real, y_pred)
varexp_train[i] <- get.rsq(Y, y.train.pred)

}

#print(paste("Found features:", ncol(id.data)))
#print(paste("Total features:", length(selected.features)-1))
#-1 because intercept term does not count

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

