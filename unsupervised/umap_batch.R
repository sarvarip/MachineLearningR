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
if (length(args) < 2) {
stop('hierarchical.R <methylation-data> <ivf-indicator>',
     call.=FALSE)
}
methylation.filename <- args[1]
ivf <- args[2]

library(umap)
library(irlba)

if (ivf == 0) {
b1 <- read.csv('../input/OriginalSamples.txt', sep='\n', header=F)
b2 <- read.csv('../input/IVFSamples.txt', sep='\n', header=F)
b3 <- read.csv('../input/AstonSamples.txt', sep='\n', header=F)

rownames(b1) <- b1$V1
rownames(b2) <- b2$V1
rownames(b3) <- b3$V1

b1$class <- 1
b2$class <- 2
b3$class <- 3

b1$V1 <- NULL
b2$V1 <- NULL
b3$V1 <- NULL

b <- rbind(b1, b2, b3)
}

df <- read.table(methylation.filename, header=T)
rownames(df) <- gsub("X", "", rownames(df))

if (ivf == 1) {
df_s <- read.csv(file="../input/ratiodat", header=TRUE, sep="\t")
df_s <- df_s[!(is.na(df_s$Female.Age) | df_s$Female.Age==""), ]
df_s <- df_s[!(is.na(df_s$Ratio) | df_s$Ratio==""), ]
rownames(df_s) <- df_s$Laboratory.ID
df_s$Laboratory.ID <- NULL
df_s$Female.Age <- NULL
data <- merge(df, df_s, by=0)
rownames(data) <- data$Row.names
data$Row.names <- NULL
data$euploid <- as.factor(as.integer((data$Ratio==0)))
data$Ratio <- NULL
train <- data
train$euploid <- NULL
train <- as.matrix(train)
} else {
df_m <- merge(df, b, by=0)
rownames(df_m) <- df_m$Row.names
df_m$Row.names <- NULL
train <- df_m
train$class <- NULL
train <- as.matrix(train)
}

set.seed(1)
df.umap <- umap(train)

filename <- sprintf('%s_UMAP.pdf', methylation.filename)
pdf(filename, width=6.5, height=8, paper='US')
if (ivf ==1){
plot(df.umap$layout, t='n', main="umap")
colors = rainbow(length(unique(data$euploid)))
names(colors) = unique(data$euploid)
text(df.umap$layout, labels=data$euploid, col=colors[data$euploid])
} else {
plot(df.umap$layout, t='n', main="umap")
colors = c('red', 'green', 'blue')
text(df.umap$layout, labels=df_m$class, col=colors[df_m$class])
}
dev.off()


if (ivf==1){
set.seed(1)
df.umap <- umap(df)
filename <- sprintf('%s_UMAP_full.pdf', methylation.filename)
pdf(filename, width=6.5, height=8, paper='US')
plot(df.umap$layout, main="umap all samples")
dev.off()
}

df.svd <- irlba::irlba(train, nv = 100)

filename <- sprintf('%s_Eigenvals.pdf', methylation.filename)
pdf(filename, width=6.5, height=8, paper='US')
plot(df.svd$d, main="Eigenvals")
dev.off()

df.reduced <- as.matrix(train) %*% df.svd$v[,1:2]

filename <- sprintf('%s_PCA.pdf', methylation.filename)
pdf(filename, width=6.5, height=8, paper='US')
plot(df.reduced, main="PCA", xlab="PC1(60%)", ylab="PC2(2%)")
dev.off()




