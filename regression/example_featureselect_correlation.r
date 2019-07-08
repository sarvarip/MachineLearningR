library(glmnet)
df <- readRDS('methylation.rds')
ages <- read.csv(file="agedat", header=TRUE, sep="\t")
rownames(df) <- gsub("X", "", rownames(df))
rownames(ages) <- ages$Laboratory.ID
ages$Laboratory.ID <- NULL
data <- merge(df, ages, by=0)
rownames(data) <- data$Row.names
data$Row.names <- NULL

rsq_list <- rep(0,5)
mse_list <- rep(0,5)
corr_list <- rep("",500)

for (i in 1:5) {

data <- data[sample(nrow(data)),]

Y_ <- data$Male.Age
X_ <- data #[,1:500]
X_$Male.Age <- NULL
cutoff = round(0.8*length(rownames(X_)))
corrs = apply(X_[1:cutoff,],2,cor,y=Y_[1:cutoff])
corrs <- corrs[order(-corrs)]
start = 1+ (i-1)*100
end = i*100
corr_list[start:end] = names(corrs[1:100])
X_red <- X_[,names(corrs[1:100])]
X <- as.matrix(X_red[1:cutoff,])
Y <- as.matrix(Y_[1:cutoff])

cv0 = cv.glmnet(X,Y,type.measure="mse",alpha=0,standardize=T)

y_real <- Y_[(cutoff+1):length(Y_)]
x_pred <- as.matrix(X_red[(cutoff+1):length(Y_),])
y_pred <- predict(cv0, newx = x_pred, s="lambda.min")
result.lm = lm(y_real ~ y_pred)
cbind(y_real, y_pred)
rsq_list[i] <- summary(result.lm)$r.squared
mse_list[i] <- mean((y_real - y_pred)^2)

}

corr_dat = as.data.frame(table(corr_list))
sorted_corr_dat = corr_dat[order(-corr_dat$Freq),]
print(sorted_corr_dat)
cbind(y_real, y_pred)
print(paste("Mean R squared: ", mean(rsq_list)))
print(paste("Std dev of R squared: ", sd(rsq_list)))
print(paste("Mean MSE: ", mean(mse_list)))
print(paste("Std dev of MSE: ", sd(mse_list)))

