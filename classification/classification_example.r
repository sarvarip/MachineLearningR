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
library(magrittr)
library(glmnet)
library(ROCR)
library(data.table)
library(Hmisc)
library(DMwR)

get_performance_vals <- function(y_real, y_pred){
    areaundercurve <- prediction(y_pred, y_real) %>% performance(measure = "auc") %>% .@y.values %>% as.numeric()
    t <- table(y_real, as.integer((y_pred > 0.5))==y_real)
    tn = t["0", "TRUE"] 
    fn = t["1", "FALSE"]
    tp = t["1", "TRUE"]
    fp = t["0", "FALSE"]

    dividecatch <- function(x,y) ifelse(y==0,0,base:::"/"(x,y))

    precision = as.numeric(dividecatch(tp,(tp+fp)))
    recall = as.numeric(dividecatch(tp,(tp+fn)))
    fscore = as.numeric(dividecatch(2*precision*recall,(precision+recall)))
    mcc = as.numeric(dividecatch((tp*tn-fp*fn),sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))))
    specificity = as.numeric(dividecatch(tn,(tn+fp)))
    accuracy = as.numeric(dividecatch((tp+tn),(tp+tn+fp+fn)))

    my_list <- list("auc" = areaundercurve, "precision" = precision, "recall" = recall, "fscore" = fscore, "mcc" = mcc, "specificity" = specificity, "accuracy" = accuracy)
    return (my_list)
}

df <- read.delim("13_methylation", header = TRUE, sep = "\t")
rownames(df) <- df$X
df$X <- NULL
df_s <- read.csv(file="filtered_successdatratio", header=TRUE, sep="\t")
df_s <- df_s[!(is.na(df_s$Female.Age) | df_s$Female.Age==""), ]
df_s <- df_s[!(is.na(df_s$Ratio) | df_s$Ratio==""), ]
rownames(df_s) <- df_s$Laboratory.ID
df_s$Laboratory.ID <- NULL
dat <- merge(df, df_s, by=0)
dat$Numsuccess <- NULL
dat$Numtrials <- NULL

f1_list <- list()
recall_list <- list()
specificity_list <- list()
auc_list <- list()
mcc_list <- list()
accuracy_list <- list()
precision_list <- list()

opt.cut = function(perf, pred){
    cut.ind = mapply(FUN=function(x, y, p){
        d = (x - 0)^2 + (y-1)^2
        ind = which(d == min(d))
        c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
            cutoff = p[[ind]])
    }, perf@x.values, perf@y.values, pred@cutoffs)
}

for (j in 1:10) {

    print("Outer iteration number: ")
    print(j)

    data <- dat[sample(nrow(dat)),]
    data$Ind = rownames(data)
    rownames(data) <- data$Row.names
    data$Row.names <- NULL
    data$euploid <- (data$Ratio==0) %>% as.integer() %>% as.factor()

    cv_fold = 5
    cv_Ind = list()
    for (i in levels(data$euploid)){
      tmp = data[data$euploid == i,c("euploid", "Ind")]
      tmp$fold <- cut2(1:dim(tmp)[1], g = cv_fold)
      levels(tmp$fold) = 1: cv_fold
      cv_Ind[[i]] = tmp
    }

    cv_Ind %<>% rbindlist() %>% setorder(Ind)
    data %>% setorder(Ind)

    data$Ratio <- NULL
    data$Ind <- NULL
    
    cv_Ind_train = cv_Ind$fold != 5
    data_train = data[cv_Ind_train,]
    cv_Ind_test = cv_Ind$fold ==5
    data_test = data[cv_Ind_test,]

    score_params = list()
    model_params = list()

    lambdas_to_try <- 10^seq(-6, -1, length.out = 10)
    for (lambda in seq(lambdas_to_try)) {

        score_fold = list()
        
        cv_fold = 5
        cv_Ind = list()
        data_train$Ind = seq(length=nrow(data_train))
        for (i in levels(data_train$euploid)){
            tmp = data_train[data_train$euploid == i,c("euploid", "Ind")]
            tmp$fold <- cut2(1:dim(tmp)[1], g = cv_fold)
            levels(tmp$fold) = 1: cv_fold
            cv_Ind[[i]] = tmp
        }
        
        cv_Ind %<>% rbindlist() %>% setorder(Ind)
        data_train %>% setorder(Ind)
        
        data_train$Ind <- NULL

        for (foldid in 1:5){

            cv_Ind_train = cv_Ind$fold != foldid
            data_cvtrain = data_train[cv_Ind_train,]
            cv_Ind_test = cv_Ind$fold == foldid
            data_cvtest = data[cv_Ind_test,]
            smote_train <- SMOTE(euploid ~ ., data  = data_train)
            #print(table(smote_train$euploid))

            Y_ <- smote_train$euploid
            X_ <- smote_train
            form <- euploid ~ .#*Female.Age
            X <- model.matrix(form, data = X_)
            Y <- Y_

            model <- glmnet(X,Y,family="binomial",alpha = 0,lambda = lambdas_to_try[lambda],standardize = TRUE)
            
            y_real <- data_cvtest$euploid
            x_pred_ <- data_cvtest
            form <- euploid ~ .#*Female.Age
            x_pred <- model.matrix(form, data = x_pred_)
            y_pred <- predict(model, newx = x_pred, type = "response")
            my_list <- get_performance_vals(y_real, y_pred)
            score_fold[foldid] = my_list$fscore

        }
        score_params[lambda] =mean(as.numeric(score_fold))
        model_params[[lambda]] = model

    }
    best_param = which.min(score_params)
    print("F1 on crossval set: ")
    print(score_params[best_param])
    chosen_lambda = lambdas_to_try[best_param]
    print("Chosen lambda: ")
    print(chosen_lambda)
    print("Index of chosen lambda: ")
    print(best_param)
    chosen_model = model_params[[best_param]]
    
    y_real <- data_test$euploid
    x_pred_ <- data_test
    form <- euploid ~ .#*Female.Age
    x_pred <- model.matrix(form, data = x_pred_)
    y_pred <- predict(chosen_model, newx = x_pred, type = "response")

    table(y_real, as.integer((y_pred > 0.5))==y_real) %>% prop.table()
    t <- table(y_real, as.integer((y_pred > 0.5))==y_real)
    print(t)

    my_list <- get_performance_vals(y_real, y_pred)
    areaundercurve = my_list$auc
    F1 = my_list$fscore
    recall = my_list$recall
    specificity = my_list$specificity
    mcc = my_list$mcc
    accuracy = my_list$accuracy
    precision = my_list$precision

    auc_list[j] <- areaundercurve
    f1_list[j] <- F1
    recall_list[j] <- recall
    specificity_list[j] <- specificity
    mcc_list[j] <- mcc
    accuracy_list[j] <- accuracy
    precision_list[j] <- precision

    print("Predictions: ")
    print(cbind(y_real, y_pred))
    print("F1 score: ")
    print(F1)
    print("Recall / Sensitivity: ")
    print(recall)
    print("Specificity: ")
    print(specificity)
    print("Area under the curve: ")
    print(areaundercurve)
    print("MCC: ")
    print(mcc)
    print("Accuracy: ")
    print(accuracy)
    print("Precision: ")
    print(precision)
    print("Best cutoff and performance values: ")
    predi <- prediction(y_pred, y_real)
    roc.perf <- prediction(y_pred, y_real) %>% performance(measure = "tpr", x.measure = "fpr")
    print(opt.cut(roc.perf, predi))

    if (j==1) {
        tmp_coeffs <- coef(chosen_model)
        df_coef <- data.frame(row.names = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
    } else {
        tmp_coeffs <- coef(chosen_model)
        df_coef2 <- data.frame(row.names = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
        df_coef <- merge(df_coef, df_coef2, by="row.names")
        rownames(df_coef) = df_coef$Row.names
        df_coef$Row.names <- NULL
        tmp <- paste("coefficient_", j)
        colnames(df_coef)[j] <- tmp
    }

}

df_coef
apply(df_coef, 1, mean)

png("AUC.png")
prediction(y_pred, y_real) %>% performance(measure = "tpr", x.measure = "fpr") %>% plot()
dev.off()

f1_list = as.numeric(f1_list)
recall_list = as.numeric(recall_list)
specificity_list= as.numeric(specificity_list)
auc_list = as.numeric(auc_list)
mcc_list = as.numeric(mcc_list)
precision_list = as.numeric(precision_list)
accuracy_list = as.numeric(accuracy_list)

print(paste("Mean F1: ", mean(f1_list), "+/-", sd(f1_list)))
print(paste("Mean recall: ", mean(recall_list), "+/-", sd(recall_list)))
print(paste("Mean specificity: ", mean(specificity_list), "+/-", sd(specificity_list)))
print(paste("Mean auc: ", mean(auc_list), "+/-", sd(auc_list)))
print(paste("Mean mcc: ", mean(mcc_list), "+/-", sd(mcc_list)))
print(paste("Mean precision: ", mean(precision_list), "+/-", sd(precision_list)))
print(paste("Mean accuracy: ", mean(accuracy_list), "+/-", sd(accuracy_list)))


