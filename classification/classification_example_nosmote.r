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

df <- read.delim("13_methylation", header = TRUE, sep = "\t")
rownames(df) <- df$X
df$X <- NULL
df_s <- read.csv(file="filtered_successdatratio", header=TRUE, sep="\t")
df_s <- df_s[!(is.na(df_s$Female.Age) | df_s$Female.Age==""), ]
df_s <- df_s[!(is.na(df_s$Ratio) | df_s$Ratio==""), ]
rownames(df_s) <- df_s$Laboratory.ID
df_s$Laboratory.ID <- NULL
dat <- merge(df, df_s, by=0)

f1_list <- rep(0,20)
recall_list <- rep(0,20)
specificity_list <- rep(0,20)
auc_list <- rep(0,20)


opt.cut = function(perf, pred){
    cut.ind = mapply(FUN=function(x, y, p){
        d = (x - 0)^2 + (y-1)^2
        ind = which(d == min(d))
        c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
            cutoff = p[[ind]])
    }, perf@x.values, perf@y.values, pred@cutoffs)
}


for (j in 1:20) {

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

    cv_Ind_train = cv_Ind$fold != 5
    cv_Ind_test = cv_Ind$fold ==5

    train_foldid = cv_Ind[cv_Ind$fold != 5]
    train_foldid = as.numeric(train_foldid$fold)

    Y_ <- data$euploid
    X_ <- data
    X_$euploid <- NULL
    X_$Ratio <- NULL
    X_$Numsuccess <- NULL
    X_$Numtrials <- NULL
    X_$Ind <- NULL
    X <- data.matrix(X_[cv_Ind_train,])
    Y <- Y_[cv_Ind_train]


    cv0 = cv.glmnet(X,Y,family = "binomial", type.measure="deviance", foldid = train_foldid, standardize=T, alpha=0)

    y_real <- Y_[cv_Ind_test]
    x_pred <- as.matrix(X_[cv_Ind_test,])
    y_pred <- predict(cv0, newx = x_pred, type = "response",  s="lambda.min")


    areaundercurve <- prediction(y_pred, y_real) %>% performance(measure = "auc") %>% .@y.values %>% as.numeric()
    table(y_real, as.integer((y_pred > 0.5))==y_real) %>% prop.table()
    t <- table(y_real, as.integer((y_pred > 0.5))==y_real)

    print(t)

    tn = t["0", "TRUE"] 
    fn = t["1", "FALSE"]
    tp = t["1", "TRUE"]
    fp = t["0", "FALSE"]

    precision = tp/(tp+fp) #1-FDR

    if (tp==0){
        recall = 0
        F1 = 0
    } else {
        recall = tp/(tp+fn) #sensitivity, TPR
        F1 = 2*precision*recall / (precision+recall)
    }
    specificity = tn/(tn+fp)

    auc_list[j] <- areaundercurve
    f1_list[j] <- F1
    recall_list[j] <- recall
    specificity_list[j] <- specificity

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
    print("Best cutoff and performance values: ")
    predi <- prediction(y_pred, y_real)
    roc.perf <- prediction(y_pred, y_real) %>% performance(measure = "tpr", x.measure = "fpr")
    print(opt.cut(roc.perf, predi))


    if (j==1) {
        tmp_coeffs <- coef(cv0, s=cv0$lambda.min)
        df_coef <- data.frame(row.names = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
    } else {
        tmp_coeffs <- coef(cv0, s=cv0$lambda.min)
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

print(paste("Mean F1: ", mean(f1_list), "+/-", sd(f1_list)))
print(paste("Mean recall: ", mean(recall_list), "+/-", sd(recall_list)))
print(paste("Mean specificity: ", mean(specificity_list), "+/-", sd(specificity_list)))
print(paste("Mean auc: ", mean(auc_list), "+/-", sd(auc_list)))

