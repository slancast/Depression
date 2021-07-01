#!/usr/bin/R
#Just what the titles says. Running c-means on the depression data

metabolomics_data <- read.table("/Users/slancast/Desktop/Projects/Depression/Metab_psych_data.csv", sep = ",", stringsAsFactors = F)
colnames(metabolomics_data) <- metabolomics_data[1,]
metabolomics_data <- metabolomics_data[-1,]

Time <- metabolomics_data$Time
metabolomics_data <- metabolomics_data[,-2] #Removing Time from the dataframe
data_table_2 <- data.frame(lapply(metabolomics_data, as.character), stringsAsFactors=FALSE, row.names = rownames(metabolomics_data)) 
data_table_2 <- data.frame(lapply(data_table_2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(metabolomics_data)) #This needs to be done for the transcript data, but causes problems with the pcl data.
data_table_2[,6:ncol(data_table_2)] <- log(data_table_2[,6:ncol(data_table_2)], 2)

data_table_3 <- aggregate(data_table_2, list(Time), FUN = mean)
rownames(data_table_3) <- data_table_3[,1]
data_table_3 <- data_table_3[,-1]
data_table_4 <- t(data_table_3 )

library(matrixStats)
library(Mfuzz)
eset <- ExpressionSet(data_table_4) #Creating the type expression set with the metadata rows as a different argument
eset <- standardise2(eset) 
expression <- exprs(eset)
eset <- ExpressionSet(na.omit(expression )) #NAs are introduced after standarizing and this is the only way to get rid of them.
eset <- standardise2(eset)
expression <- exprs(eset)
m = mestimate(eset)


cluster_number <- Dmin(eset,m,crange=seq(2,15,1),repeats=20,visu=TRUE)
mfuzzcl <- mfuzz(eset, c=15, m=m)
library(MASS)
cluster_labels <- mfuzzcl$cluster

xaxis_ticks = c("T1","T2","T3","T4","T5")

pdf(paste("~/Desktop/metabolomics_depression.pdf",sep=""))
par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
    mar = c(2,1,1,1) + 0.1)
mfuzz.plot2(eset,cl=mfuzzcl,ylim=c(-3,3),mfrow=c(3,3), time.labels = xaxis_ticks, ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.main=1.1,cex.lab=0.1)
dev.off()

write.table(cluster_labels, file = "~/Desktop/metabolomics_depression_c-means.txt", sep = "\t", quote = F)

