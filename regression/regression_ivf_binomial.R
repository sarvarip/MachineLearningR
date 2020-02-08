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

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
stop('regression.R <methylation-data> <ivf-data> <alpha>',
     call.=FALSE)
}
methylation.filename <- args[1]
metadata.filename <- args[2]
alpha.lasso <- args[3]

library(glmnet)
df <- read.table(methylation.filename, header=T)
rownames(df) <- gsub("X", "", rownames(df))
if (length(grep(".rds", metadata.filename))>0){
df_s <- readRDS(metadata.filename)
} else {
df_s <- read.csv(file=metadata.filename, header=TRUE, sep="\t", stringsAsFactors=FALSE)
df_s <- df_s[!(is.na(df_s$Female.Age) | df_s$Female.Age==""), ]
df_s$Female.Age <- as.numeric(df_s$Female.Age)
df_s <- df_s[!(is.na(df_s$succ) | df_s$succ==""), ]
df_s <- df_s[!(is.na(df_s$nosucc) | df_s$nosucc==""), ]
}
if ("Laboratory.ID" %in% colnames(df_s)) {
rownames(df_s) <- df_s$Laboratory.ID
df_s$Laboratory.ID <- NULL
}

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
print("Number of samples before row expansion")
print(n.individuals)
parts <- sample(rep(1:cross.val.fold, length.out=n.individuals))

accuracy_list <- rep(0,cross.val.fold)
base_accuracy <- rep(0, cross.val.fold)
accuracy_train <- rep(0,cross.val.fold)
accuracy_majority <- rep(0, cross.val.fold)
#Theoretical maximum accuracy
accuracy_max <- rep(0, cross.val.fold)
varexp_train <- rep(0,cross.val.fold)
varexp_list <- rep(0,cross.val.fold)
if (alpha.lasso > 0) {
zeroed_list <- rep(0,cross.val.fold)
}

for (i in 1:cross.val.fold) {

training.individuals <- rownames(data)[which(parts != i)]
testing.individuals <- rownames(data)[which(parts == i)]
X_ <- data[training.individuals,]

#Expanding rows
X_$trials <- X_$succ + X_$nosucc
X_.expanded <- X_[rep(row.names(X_), X_$trials), ]
succdat <- X_[, c("nosucc", "succ")]
success <- rep(rep(c(0,1), nrow(succdat)), c(t(succdat)))
X_.expanded$outcome <- success
print(paste("Training: Class 1 has ", sum(success==1), "number of samples"))
print(paste("Training: Class 0 has ", sum(success==0), "number of samples"))
X_.expanded$succ <- NULL
X_.expanded$nosucc <- NULL
X_.expanded$trials <- NULL
X_ <- X_.expanded

Y <- data.matrix(X_$outcome)
base <- X_[, c("Female.Age", "outcome")]
X_$outcome <- NULL
X <- data.matrix(X_)

cv0 = cv.glmnet(X,Y,type.measure="deviance",alpha=alpha.lasso,
                standardize=T, family="binomial")
if (alpha.lasso > 0) {
cs <- coef(cv0, s='lambda.min')
zeroed <- sum(cs == 0)
zeroed_list[i] <- zeroed
}
basefit <- glm(outcome ~ Female.Age, data=base, family="binomial")

#Calculating theoretically maximum accuracy
X.test <- data[testing.individuals,]
check <- X.test[, c("nosucc", "succ")]
check$max <- apply(check, 1, max)
check$total <- check$succ + check$nosucc
accuracy_max[i] <- sum(check$max) / sum(check$total)

#Expanding rows
X.test$trials <- X.test$succ + X.test$nosucc
X.test.expanded <- X.test[rep(row.names(X.test), X.test$trials), ]
succdat <- X.test[, c("nosucc", "succ")]
success <- rep(rep(c(0,1), nrow(succdat)), c(t(succdat)))
X.test.expanded$outcome <- success
number_succ <- sum(success==1)
number_nosucc <- sum(success==0)
accuracy_majority[i] <- max(number_succ, number_nosucc) / (number_succ+number_nosucc)
print(paste("Testing: Class 1 has ", number_succ, "number of samples"))
print(paste("Testing: Class 0 has ", number_nosucc, "number of samples"))
X.test.expanded$succ <- NULL
X.test.expanded$nosucc <- NULL
X.test.expanded$trials <- NULL
X.test <- X.test.expanded


y_real <- data.matrix(X.test$outcome)
base.test <- X.test[, c("Female.Age", "outcome")]
X.test$outcome <- NULL
x_pred <- data.matrix(X.test)

y_pred <- predict(cv0, newx = x_pred, type = "response", s="lambda.min")
y.train.pred <- predict(cv0, newx = X, type = "response", s="lambda.min")
y.baseline.pred <- predict(basefit, newdata = base.test, type= "response")
print(cbind(y_real, y_pred))
accuracy_list[i] <- mean(as.integer((y_pred > 0.5))==y_real)
base_accuracy[i] <- mean(as.integer((y.baseline.pred > 0.5))==y_real)
print(paste("Testing accuracy in this iteration: ", accuracy_list[i]))
print(paste("Majority prediction testing accuracy: ", accuracy_majority[i]))
print(paste("Theoretical maximum testing accuracy: ", accuracy_max[i]))
print(paste("Female age only accuracy: ", base_accuracy[i]))
accuracy_train[i] <- mean(as.integer((y.train.pred > 0.5))==Y)
print(paste("Training accuracy in this iteration: ", accuracy_train[i]))
varexp_list[i] <- get.rsq(y_real, y_pred)
varexp_train[i] <- get.rsq(Y, y.train.pred)

}

print(paste("Mean accuracy: ", mean(accuracy_list)))
print(paste("Std dev of accuracy: ", sd(accuracy_list)))
print(paste("Mean majority prediction accuracy: ", mean(accuracy_majority)))
print(paste("Std dev of majority accuracy: ", sd(accuracy_majority)))
print(paste("Mean theoretical maximum accuracy: ", mean(accuracy_max)))
print(paste("Std dev of maximum accuracy: ", sd(accuracy_max)))
print(paste("Mean female age accuracy: ", mean(base_accuracy)))
print(paste("Std dev of female age accuracy: ", sd(base_accuracy)))
print(paste("Training set: Mean accuracy: ", mean(accuracy_train)))
print(paste("Training set: Std dev of accuracy: ", sd(accuracy_train)))
print(paste("Mean of variance explained: ", mean(varexp_list)))
print(paste("Std dev of variance explained: ", sd(varexp_list)))
print(paste("Training set: Mean of variance explained: ", 
            mean(varexp_train)))
print(paste("Training set: Std dev of variance explained: ", 
            sd(varexp_train)))
print(paste("Number of features: ", ncol(data)-1))
if (alpha.lasso > 0) {
print(paste("Average number of zeroed features: ",
            mean(zeroed_list)))
}



