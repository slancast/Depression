#!/usr/bin/R

#Reading in the data table generated by DADA2 from the program data2.R 
dada2_data_table <- read.table("/Users/slancast/Box/School for The Work Data -- Sam/Microbiome/depression_dada2.txt", stringsAsFactors = F)
duplicated_sample1 <- which(dada2_data_table[,1] %in% "559299082_NA0021494110_5_gut") #This sample did not work and is duplicated in the subsequent row, so I'm taking it out
dada2_data_table <- dada2_data_table[-duplicated_sample1,]

#Finding the relative abundances
#getting the numeric values out of the data table
dada2_data_table2 <- dada2_data_table2[8:nrow(dada2_data_table2),2:ncol(dada2_data_table2)]
dada2_data_table2 <- data.frame(lapply(dada2_data_table, as.character), stringsAsFactors=FALSE, row.names = rownames(dada2_data_table)) 
dada2_data_table2 <- data.frame(lapply(dada2_data_table, as.numeric), stringsAsFactors=FALSE, row.names = rownames(dada2_data_table))
row_summation <- rowSums(dada2_data_table2)
for (i in 1:ncol(dada2_data_table2)){
  dada2_data_table2[,i] <- dada2_data_table2[,i]/row_summation
}
dada2_data_table[8:nrow(dada2_data_table),2:ncol(dada2_data_table)] <- dada2_data_table2
write.table(dada2_data_table, file = paste0("/Users/slancast/Desktop/depression_dada2_relative.txt"),sep="\t", quote=F, col.names = F)

