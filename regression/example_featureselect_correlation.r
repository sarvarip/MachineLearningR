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

df <- readRDS('../methylation.rds')
rownames(df) <- gsub("X", "", rownames(df))

df_s <- read.csv(file="../agedat", header=TRUE, sep="\t")
df_s <- df_s[!(is.na(df_s$Male.Age) | df_s$Male.Age==""), ]
rownames(df_s) <- df_s$Laboratory.ID
df_s$Laboratory.ID <- NULL

data <- merge(df, df_s, by=0)
rownames(data) <- data$Row.names
data$Row.names <- NULL

data$Age <- data$Male.Age
data$Male.Age <- NULL

get.rsq <- function(y, y.hat) {
  sst <- sum((y - mean(y))^2)
  sse <- sum((y.hat - y)^2)
  return(1 - sse/sst)
}

cross.val.fold <- 5
n.individuals <- nrow(data)
print("Number of samples")
print(n.individuals)
parts <- sample(rep(1:cross.val.fold, length.out=n.individuals))
rsq_list <- rep(0,cross.val.fold)
varexp_list <- rep(0,cross.val.fold)
mse_list <- rep(0,cross.val.fold)
num_features <- 10
hyperparams.to.try <- 5
corr_list <- rep("", cross.val.fold*cross.val.fold*num_features)
data$Numsuccess <- NULL
data$Numtrials <- NULL

for (j in 1:cross.val.fold) {

data <- data[sample(nrow(data)),]

training.individuals <- rownames(data)[which(parts != j)]
testing.individuals <- rownames(data)[which(parts == j)]
data_train <- data[training.individuals,]
data_test <- data[testing.individuals,]

score_params <- rep(0,cross.val.fold)
model_params = list()
corr_paramslist = rep("", cross.val.fold*hyperparams.to.try*num_features)

lambdas_to_try <- 10^seq(-1, 4, length.out = hyperparams.to.try)
for (lambda in seq(lambdas_to_try)) {
    score_fold <- rep(0,cross.val.fold)

    n.individuals <- nrow(data_train)
    parts <- sample(rep(1:cross.val.fold, length.out=n.individuals))

    for (i in 1:cross.val.fold) {
        print("Info")
        print(lambda)
        print(i)
        
        training.individuals <- rownames(data_train)[which(parts != i)]
        testing.individuals <- rownames(data_train)[which(parts == i)]
        data_train_train <- data_train[training.individuals,]
        X_ <- data_train_train
        X_$Ratio <- NULL

        corrs = apply(X_,2,cor,y=data_train_train$Ratio)
        corrs <- corrs[order(-corrs)]
        start = 1+(lambda-1)*cross.val.fold*num_features+(i-1)*num_features
        end = (lambda-1)*cross.val.fold*num_features+i*num_features
        corr_paramslist[start:end] = names(corrs[1:num_features])

        X_red <- X_[,names(corrs[1:num_features])]

        X <- as.matrix(X_red)
        Y <- as.matrix(data_train_train$Ratio)
        
        model <- glmnet(X,Y,alpha = 0,lambda = lambdas_to_try[lambda],standardize = TRUE)
        data_train_test <- data_train[testing.individuals,]
        y_real <- as.matrix(data_train_test$Ratio)
        x_pred <- data_train_test[,names(corrs[1:num_features])]
        x_pred$Ratio <- NULL
        x_pred <- as.matrix(x_pred)
        y_pred <- predict(model, newx = x_pred)
        score_fold[i] <- get.rsq(y_real, y_pred)
    }
    score_params[lambda] = mean(score_fold)
    model_params[[lambda]] = model
}
best_param = which.max(score_params)
print("Variance explained on crossval set: ")
print(score_params[best_param])
chosen_lambda = lambdas_to_try[best_param]
print("Chosen lambda: ")
print(chosen_lambda)
print("Index of chosen lambda: ")
print(best_param)
chosen_model = model_params[[best_param]]
start = 1+ (j-1)*cross.val.fold*num_features
end = j*cross.val.fold*num_features
start_ = 1+(best_param-1)*cross.val.fold*num_features
end_ = best_param*cross.val.fold*num_features
corr_list[start:end] = corr_paramslist[start_:end_]
corr_dat = as.data.frame(table(corr_list))
corr_dat = corr_dat[!(corr_dat$corr_list==""),]
sorted_corr_dat = corr_dat[order(-corr_dat$Freq),]
features <- sorted_corr_dat[1:num_features,"corr_list"]

y_real <- as.numeric(data_test$Ratio)
x_pred <- data_test #[,1:500]
x_pred <- x_pred[,features]
x_pred$Ratio <- NULL
x_pred <- as.matrix(x_pred)

y_pred <- predict(chosen_model, newx = x_pred)
result.lm = lm(y_real ~ y_pred)
cbind(y_real, y_pred)
rsq_list[j] <- summary(result.lm)$r.squared
mse_list[j] <- mean((y_real - y_pred)^2)
varexp_list[j] <- get.rsq(y_real, y_pred)

}

corr_dat = as.data.frame(table(corr_list))
sorted_corr_dat = corr_dat[order(-corr_dat$Freq),]
features <- sorted_corr_dat[1:num_features,]
saveRDS(features, 'features_10_corrselect_ivfdat.rds') 
# in the future will also need to save lambda
cbind(y_real, y_pred)
print(paste("Mean R squared: ", mean(rsq_list)))
print(paste("Std dev of R squared: ", sd(rsq_list)))
print(paste("Mean MSE: ", mean(mse_list)))
print(paste("Std dev of MSE: ", sd(mse_list)))
print(paste("Mean of variance explained: ", mean(varexp_list)))
print(paste("Std dev of variance explained: ", sd(varexp_list)))

