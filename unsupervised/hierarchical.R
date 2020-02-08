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
if (length(args) < 1) {
stop('hierarchical.R <methylation-data>',
     call.=FALSE)
}
methylation.filename <- args[1]


library(DMwR2)

df <- read.table(methylation.filename, header=T)
rownames(df) <- gsub("X", "", rownames(df))
o <- outliers.ranking(df)
print(cbind(o$rank.outliers, o$prob.outliers[o$rank.outliers]))
hc <- hclust(dist(df), method="complete") #no scaling, since 
#methylation levels have same units

filename <- sprintf('%s_hclust.pdf', methylation.filename)
pdf(filename, width=6.5, height=8, paper='US')
plot(hc)
dev.off()





