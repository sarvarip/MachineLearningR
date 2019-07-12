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

df <- readRDS("age_1000topvariance.rds")
#rownames(df) <- gsub("X", "", rownames(df))
df_s <- read.csv(file="ages_successdatratio", header=TRUE, sep="\t")
df_s <- df_s[!(is.na(df_s$Female.Age) | df_s$Female.Age==""), ]
df_s <- df_s[!(is.na(df_s$Ratio) | df_s$Ratio==""), ]

data <- merge(df, df_s, by=0)
rownames(data) <- data$Row.names
data$Row.names <- NULL

rsq_list <- rep(0,5)
mse_list <- rep(0,5)
num_features = 10
corr_list <- rep("", 25*num_features)

for (j in 1:5) {

data <- data[sample(nrow(data)),]

data$Numsuccess <- NULL
data$Numtrials <- NULL
cutoff = round(0.8*nrow(data))
data_train <- data[1:cutoff,]
data_test <- data[(cutoff+1):nrow(data),]

score_params <- rep(0,5)
model_params = list()
corrs_params = list()
corr_paramslist = rep("", 50*num_features)
cutoff_outer <- as.integer(cutoff)

lambdas_to_try <- 10^seq(-6, -1, length.out = 10)
for (lambda in seq(lambdas_to_try)) {
    score_fold <- rep(0,5)
    for (i in 1:5) {
        print("Info")
        print(lambda)
        print(i)
        Y_ <- data_train$Ratio
        X_ <- data_train #[,1:500]
        X_$Ratio <- NULL
        cutoff <- as.integer(round(0.8*nrow(data_train)))
        corrs = apply(X_[1:cutoff,],2,cor,y=Y_[1:cutoff])
        corrs <- corrs[order(-corrs)]
        start = 1+(lambda-1)*5*num_features+(i-1)*num_features
        end = (lambda-1)*5*num_features+i*num_features
        corr_paramslist[start:end] = names(corrs[1:num_features])
        X_red <- X_[,names(corrs[1:num_features])]
        X <- as.matrix(X_red[1:cutoff,])
        Y <- as.matrix(Y_[1:cutoff])
        model <- glmnet(X,Y,alpha = 0,lambda = lambdas_to_try[lambda],standardize = TRUE)
        y_real <- Y_[(cutoff+1):cutoff_outer]
        x_pred <- X_red[(cutoff+1):cutoff_outer,]
        x_pred <- as.matrix(x_pred)
        y_pred <- predict(model, newx = x_pred)
        result.lm = lm(y_real ~ y_pred)
        score_fold[i] <- summary(result.lm)$r.squared
    }
    score_params[lambda] = mean(score_fold)
    model_params[[lambda]] = model
    corrs_params[[lambda]] = names(corrs[1:num_features])
}
best_param = which.min(score_params)
print("R^2 on crossval set: ")
print(score_params[best_param])
chosen_lambda = lambdas_to_try[best_param]
print("Chosen lambda: ")
print(chosen_lambda)
print("Index of chosen lambda: ")
print(best_param)
chosen_model = model_params[[best_param]]
features = corrs_params[[best_param]]
start = 1+ (j-1)*5*num_features
end = j*5*num_features
start_ = 1+(best_param-1)*5*num_features
end_ = best_param*5*num_features
corr_list[start:end] = corr_paramslist[start_:end_]
print(corr_list)

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

}

corr_dat = as.data.frame(table(corr_list))
sorted_corr_dat = corr_dat[order(-corr_dat$Freq),]
print(sorted_corr_dat)
cbind(y_real, y_pred)
print(paste("Mean R squared: ", mean(rsq_list)))
print(paste("Std dev of R squared: ", sd(rsq_list)))
print(paste("Mean MSE: ", mean(mse_list)))
print(paste("Std dev of MSE: ", sd(mse_list)))

