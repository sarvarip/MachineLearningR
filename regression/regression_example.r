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

rsq_list <- rep(0,20)
mse_list <- rep(0,20)

for (i in 1:20) {

data <- data[sample(nrow(data)),]
Y_ <- data$Ratio
X_ <- data
X_$Ratio <- NULL
X_$Numsuccess <- NULL
X_$Numtrials <- NULL
cutoff = round(0.8*length(rownames(X_)))
X <- as.matrix(X_[1:cutoff,])
Y <- as.matrix(Y_[1:cutoff])

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

y_real <- Y_[(cutoff+1):length(Y_)]
x_pred <- as.matrix(X_[(cutoff+1):length(Y_),])
y_pred <- predict(cv0, newx = x_pred, s="lambda.min")
result.lm = lm(y_real ~ y_pred)
cbind(y_real, y_pred)
rsq_list[i] <- summary(result.lm)$r.squared
mse_list[i] <- mean((y_real - y_pred)^2)

}

df_coef
apply(df_coef, 1, mean)
print(paste("Mean R squared: ", mean(rsq_list)))
print(paste("Std dev of R squared: ", sd(rsq_list)))
print(paste("Mean MSE: ", mean(mse_list)))
print(paste("Std dev of MSE: ", sd(mse_list)))


