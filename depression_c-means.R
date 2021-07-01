#!/usr/bin/R

dada2_data_table <- read.table("/Users/slancast/Desktop/Projects/Depression/Data/depression_dada2_relative.txt", sep="\t", stringsAsFactors = F, row.names = 1)
#Pulling out samplenames to identify kit ids
sample_name <- dada2_data_table[,1]
sample_name <- strsplit(as.character(sample_name), "_") 
samples <- c()
#Pulling out the kid IDs from the sample naems
for (i in sample_name) {
  samples <- c(samples, i[1])
}
#finding duplicated samples and averaging them
duplicated_samples <- samples[duplicated(samples)]
#Creating a blank data frame
replacement_samples <- dada2_data_table[FALSE,] #Creates an empty data frame with no rows
for (i in duplicated_samples){
  duplicated_sample <- which(samples %in% i)
  new_row <- (as.numeric(dada2_data_table[duplicated_sample[1],2:ncol(dada2_data_table)]) + as.numeric(dada2_data_table[duplicated_sample[2],2:ncol(dada2_data_table)])) / 2 #Finding the average of the duplicated rows. Obviously only having two duplicated rows is hard coded in here
  new_row <- c(i, new_row)
  names(new_row) <- colnames(dada2_data_table)
  new_row <- data.frame(t(new_row))
  replacement_samples <- rbind(replacement_samples, new_row)
}
#removing duplciates
duplicated_sample_numbers <- which(samples %in% duplicated_samples)
dada2_data_table <- dada2_data_table[-duplicated_sample_numbers,]
dada2_data_table <- rbind(dada2_data_table,replacement_samples)

dada2_metadata <- read.table("/Users/slancast/Box/School for The Work Data -- Sam/Microbiome/metadata_full_deduplicated.csv", sep=",", header=T, stringsAsFactors = F)

sample_name <- dada2_data_table[,1]
sample_name <- strsplit(as.character(sample_name), "_") 
samples <- c()
for (i in sample_name) {
  samples <- c(samples, i[1])
}
samples_order <- order(samples)
reordered_samples <- samples[samples_order]
#Remaking samples with the new data frame
samples_sequenced <- which(reordered_samples %in% dada2_metadata$kitid)

#This will only keep the rows of the samples that are sequenced.
#But it also 
dada2_data_table_reodered <- dada2_data_table[samples_order,]
data_table_2 <- dada2_data_table_reodered[samples_sequenced,]
rownames(data_table_2) <- data_table_2[,1]
data_table_2 <- data_table_2[,-1]
data_table_2 <- data.frame(lapply(data_table_2, as.character), stringsAsFactors=FALSE, row.names = rownames(data_table_2)) 
data_table_2 <- data.frame(lapply(data_table_2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(data_table_2)) #This needs to be done for the transcript data, but causes problems with the pcl data.


metadata_order <- order(dada2_metadata$kitid)
reordered_metadata <- dada2_metadata[metadata_order,]
#reversing it to save the metadata that are in samples
metadata_needed <- which(reordered_metadata$kitid %in% samples)
metadata_2 <- reordered_metadata[metadata_needed,]

data_table_3 <- aggregate(data_table_2, list(metadata_2$Time), FUN = mean)
rownames(data_table_3) <- data_table_3[,1]
data_table_3 <- data_table_3[,-1]
data_table_4 <- t(data_table_3 )

names <- data.frame(t(dada2_data_table[1:7,]))
names2 <- paste(names[,1], names[,2], names[,3], names[,4], names[,5], names[,6], names[,7],sep = "_")
names2 <- names2[-1]
names2 <- make.names(names2, unique = T)
rownames(data_table_4) <- names2

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
mfuzzcl <- mfuzz(eset, c=10, m=m)
library(MASS)
cluster_labels <- mfuzzcl$cluster

xaxis_ticks = c("T1","T3","T4","T5")


# pdf(paste("~/Desktop/depression.pdf",sep=""))
# par(mfrow=c(5,4),oma = c(5,4,1,0) + 0.1,
#     mar = c(2,1,1,1) + 0.1)
# mfuzz.plot2(eset,cl=mfuzzcl,ylim=c(-3,3),mfrow=c(3,3), time.labels = xaxis_ticks, ax.col="red",col.lab="black",col.main="green",col.sub="blue",col="blue", Xwidth=20, Xheight=20,colo="fancy", lwd=1,ylab='',xlab='',x11=FALSE,cex.main=1.1,cex.lab=0.1)
# dev.off()

# write.table(cluster_labels, file = "~/Desktop/depression_c-means.txt", sep = "\t", quote = F)



colnames(data_table_2) <- names2 
correlation_mtarix_cols <- which(colSums(data_table_2) > 0.001) #Any more than this and the pvalue hist looks weird
correlation_mtarix <- data_table_2[,correlation_mtarix_cols]

metabolomics_data <- read.table("/Users/slancast/Desktop/Projects/Depression/Metab_psych_data.csv", sep = ",", stringsAsFactors = F)
colnames(metabolomics_data) <- metabolomics_data[1,]
metabolomics_data <- metabolomics_data[-1,]

psych_data <- read.table("/Users/slancast/Desktop/Psychometric Data.csv", sep = ",", stringsAsFactors = F)
colnames(psych_data) <- psych_data[1,]
psych_data <- psych_data[-1,]

#Matching up the time and ID for the microbiome metadata to match with the other data
microbiome_time_by_id <- as.character(interaction(metadata_2$id, metadata_2$Time))
metabolome_time_by_id <- as.character(interaction(metabolomics_data$id, metabolomics_data$Time))
psych_time_by_id <- as.character(interaction(psych_data$id, psych_data$Time))

#Making these the rownames
data_table_3 <- correlation_mtarix[-which(is.na(microbiome_time_by_id )),]
microbiome_time_by_id <- na.omit(microbiome_time_by_id)
rownames(data_table_3) <- microbiome_time_by_id 
rownames(metabolomics_data) <- metabolome_time_by_id 
rownames(psych_data) <- psych_time_by_id
psych_data$depressed <- as.numeric(as.factor(psych_data$depressed))

#Transposing the datasets
data_table_4 <- data.frame(t(data_table_3))
data_table_4  <- data.frame(lapply(data_table_4 , as.character), stringsAsFactors=FALSE, row.names = rownames(data_table_4 )) #This needs to be done for the transcript data, but causes problems with the pcl data.
data_table_4  <- data.frame(lapply(data_table_4 , as.numeric), stringsAsFactors=FALSE, row.names = rownames(data_table_4 )) #This needs to be done for the transcript data, but causes problems with the pcl data.
data_table_4  <- data.frame(na.omit(data_table_4))
#data_table_4  <- data_table_4[-rows_microbes_with_NAs,] #Getting rid of the microbes with NAs
metabolomics_data_2 <- data.frame(t(metabolomics_data[,])) #-c(1:6) this will get rid of the psychometric data columns from metabolomics data if inserted
metabolomics_data_2  <- data.frame(lapply(metabolomics_data_2 , as.character), stringsAsFactors=FALSE, row.names = rownames(metabolomics_data_2)) #This needs to be done for the transcript data, but causes problems with the pcl data.
metabolomics_data_2  <- data.frame(lapply(metabolomics_data_2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(metabolomics_data_2 )) #This needs to be done for the transcript data, but causes problems with the pcl data.
metabolomics_data_2 <- data.frame(na.omit(metabolomics_data_2))
metabolomics_data_2 <- metabolomics_data_2[-c(1:5),]
metabolomics_data_2 <- log(metabolomics_data_2,2)
psych_data_2 <- data.frame(t(psych_data))
psych_data_2 <- data.frame(lapply(psych_data_2, as.character), stringsAsFactors=FALSE, row.names = rownames(psych_data_2)) #This needs to be done for the transcript data, but causes problems with the pcl data.
psych_data_2 <- data.frame(lapply(psych_data_2, as.numeric), stringsAsFactors=FALSE, row.names = rownames(psych_data_2)) #This needs to be done for the transcript data, but causes problems with the pcl data.
psych_data_2 <- psych_data_2[, which(colMeans(!is.na(psych_data_2)) > 0.8)]
library(impute) 
psych_toimpute.i <- impute.knn(as.matrix(psych_data_2)) #Imputing the individual entries missing 
psych_imputed <- data.frame(psych_toimpute.i$data) #pulling out the data from the imputed object
psych_data_2 <- data.frame(na.omit(psych_imputed))

#Here I will compile a list of metabolites that were specifically correlated with depression scores
#This was done by Daisy
#Using this targeted list of metabolites will increase statistical power, by reducing FDR penalties
Depression_Negative_MetabolomicsAnnotated <- read.table("/Users/slancast/Box/Tibshirani Lab/Metabolites/NEWDATA_DrugOut/NEWDATA_DrugOut/Annotated/Depression_Negative_MetabolomicsAnnotated.csv", sep=",", header=T, stringsAsFactors = F)
Depression_Positive_MetabolomicsAnnotated <- read.table("/Users/slancast/Box/Tibshirani Lab/Metabolites/NEWDATA_DrugOut/NEWDATA_DrugOut/Annotated/Depression_Positive_MetabolomicsAnnotated.csv", sep=",", header=T, stringsAsFactors = F)
Safety_Negative_MetabolomicsAnnotated <- read.table("/Users/slancast/Box/Tibshirani Lab/Metabolites/NEWDATA_DrugOut/NEWDATA_DrugOut/Annotated/Depression_Positive_MetabolomicsAnnotated.csv", sep=",", header=T, stringsAsFactors = F)
Safety_Positive_MetabolomicsAnnotated <- read.table("/Users/slancast/Box/Tibshirani Lab/Metabolites/NEWDATA_DrugOut/NEWDATA_DrugOut/Annotated/Depression_Positive_MetabolomicsAnnotated.csv", sep=",", header=T, stringsAsFactors = F)
# metabolite_depression_association_without_on_drug <- read.table("/Users/slancast/Box/Tibshirani Lab/Metabolites/NEWDATA_DrugOut/NEWDATA_DrugOut/results_0330_metabolites_without_drug/metabolite_depression_association_without_on_drug.txt", sep=" ", header=F, stringsAsFactors = F)
# metabolite_depression_slope_without_on_drug <- read.table("/Users/slancast/Box/Tibshirani Lab/Metabolites/NEWDATA_DrugOut/NEWDATA_DrugOut/results_0330_metabolites_without_drug/metabolite_depression_slope_without_on_drug.txt", sep=" ", header=F, stringsAsFactors = F)
# metabolite_safety_association_without_on_drug <- read.table("/Users/slancast/Box/Tibshirani Lab/Metabolites/NEWDATA_DrugOut/NEWDATA_DrugOut/results_0330_metabolites_without_drug/metabolite_safety_association_without_on_drug.txt", sep=" ", header=F, stringsAsFactors = F)
# metabolite_safety_slope_without_on_drug <- read.table("/Users/slancast/Box/Tibshirani Lab/Metabolites/NEWDATA_DrugOut/NEWDATA_DrugOut/results_0330_metabolites_without_drug/metabolite_safety_slope_without_on_drug.txt", sep=" ", header=F, stringsAsFactors = F)

#Pulling out the interesting metabolites
DNMA <- Depression_Negative_MetabolomicsAnnotated$Variables
DPMA <- Depression_Positive_MetabolomicsAnnotated$Variables
SNMA <- Safety_Negative_MetabolomicsAnnotated$Variables
SPMA <- Safety_Positive_MetabolomicsAnnotated$Variables

#Combining into single list
metabolites_of_interest <- c(DNMA, DPMA, SNMA, SPMA)
metabolites_of_interest <- unique(metabolites_of_interest )
metabolites_of_interest <- gsub("m\\.z", "m\\/z", metabolites_of_interest)

#Various other metabolites of interest. Quick correlation with the following measurements
# M1.75_237.1492m/z #gad7
# M8.57_307.1186m/z, M7.58_189.0692m/z, M14.24_232.1403m/z, M8.81_329.2332m/z, M10.60_205.1181m/z #bdi9
metabolites_of_interest2 <- c("M1.75_237.1492m/z", "M8.57_307.1186m/z", "M7.58_189.0692m/z", "M14.24_232.1403m/z", "M8.81_329.2332m/z", "M10.60_205.1181m/z")
metabolite_psych_data <- colnames(metabolomics_data)[1:6]

# Subsetting metabolites. Toggle this on and off to subset or not
total_metabolites_of_interest <- c(metabolites_of_interest, metabolites_of_interest2, metabolite_psych_data)
total_metabolites_of_interest_rows <- which(rownames(metabolomics_data_2) %in% total_metabolites_of_interest)
#metabolomics_data_2 <- metabolomics_data_2[total_metabolites_of_interest_rows,]


#Combining the datasets
library(plyr)
combined_df <- rbind.fill(psych_data_2, data_table_4) #metabolomics_data_2, 
rownames(combined_df) <- c(rownames(psych_data_2), rownames(data_table_4)) #rownames(metabolomics_data_2), 

#removing_NAs
combined_df <- t(combined_df )
combined_df2 <- na.omit(combined_df )



#Sample correlation network code
library("Hmisc")
cor <- rcorr(format(combined_df,digits=20), type="spearman")
cor.data <- cor$r
rownames(cor.data) <- colnames(combined_df)
cor.data[upper.tri(cor.data, diag = T)] <- 0
pval.data <- cor$P
pval.data[upper.tri(pval.data, diag = T)] <- NA
FDR.data <- apply(pval.data,2,p.adjust,method="BH", n = length(pval.data))
pdf(paste("~/Desktop/depression_pval_hist.pdf",sep=""))
hist(pval.data, breaks = 100, col="darkblue")
dev.off()
pdf(paste("~/Desktop/depression_cor_hist.pdf",sep=""))
hist(cor.data, breaks = 10, col="red")
dev.off()
cor.data[FDR.data > 0.05]=0
write.table(cor.data,file=paste("~/Desktop/depression_cor_data.txt",sep=""), sep="\t")
write.table(FDR.data,file=paste("~/Desktop/depression_fdr_data.txt",sep=""), sep="\t")

library("igraph")
network=graph.adjacency(abs(cor.data), weighted=T, mode="undirected", diag=F)

#The following lines will find if  the correlations are negative or positive
#and then add the properly colored edge
correlation_values <- c()
edges <- data.frame(get.edgelist(network))
for (i in 1:nrow(edges)){
  j <- as.matrix(edges[i,])
  correlation_value <- cor.data[j[2],j[1]]
  correlation_values <- c(correlation_values,correlation_value)
}
E(network)$color <- ifelse(correlation_values < 0, "red","darkblue")

V(network)$color <- c(rep("lightblue",nrow(psych_data_2)),rep("green",nrow(metabolomics_data_2)), rep("blue",nrow(data_table_4))) #


ly <- layout_with_fr(network,dim=2,grid="nogrid",niter=1000)
#ly <- layout.fruchterman.reingold(subgraph,dim=2,grid="nogrid")

pdf("/Users/slancast/Desktop/depression_network.pdf")
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(network,
     vertex.size=2.5,
     #vertex.color=V(subgraph)$color,
     vertex.label.cex=0.25,
     vertex.label.color="black",
     vertex.frame.color="black",
     vertex.label = NA,
     layout = ly,
     edge.width=abs(E(network)$weight)
)
dev.off()


if (FALSE) {
# Finding the rows that can't create correlations
# Haven't flushed out why, and might want to 
# For sure want to eliinate them though
row_to_toss1 <- (rowSums(is.na(cor.data[,1:2])) > 0)
row_to_toss2 <- which(isTRUE(row_to_toss))

microbes_with_NAs <- c("Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Muribaculaceae_NA_NA",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Muribaculaceae_NA_NA.3",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Prevotellaceae_NK3B31_group_NA.3",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Prevotellaceae_UCG.001_NA.4",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Barnesiellaceae_Barnesiella_NA.3",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Paraprevotella_NA.4",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Prevotella_2_NA.1",
"Bacteria_Firmicutes_Negativicutes_Selenomonadales_Acidaminococcaceae_Acidaminococcus_NA",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Faecalibacterium_NA.12",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Prevotellaceae_UCG.001_NA.5",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Faecalibacterium_NA.15",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Ruminococcus_1_NA.9",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Muribaculaceae_NA_NA.13",
"Bacteria_Firmicutes_Negativicutes_Selenomonadales_Veillonellaceae_Megasphaera_NA",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_NA_NA.20",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Peptostreptococcaceae_Peptoclostridium_NA",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Muribaculaceae_Muribaculum_NA",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Alloprevotella_NA.2",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Faecalibacterium_NA.17",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Barnesiellaceae_Barnesiella_NA.7",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_Butyrivibrio_NA.1",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Muribaculaceae_NA_NA.16",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_NA_NA.29",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_NA_NA.30",
"Bacteria_Fusobacteria_Fusobacteriia_Fusobacteriales_Fusobacteriaceae_Fusobacterium_NA.4",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Bacteroidaceae_Bacteroides_NA.39",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_Butyrivibrio_NA.2",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Alloprevotella_NA.3",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Prevotellaceae_NK3B31_group_NA.5",
"Bacteria_Actinobacteria_Coriobacteriia_Coriobacteriales_Atopobiaceae_Coriobacteriaceae_UCG.003_NA",
"Bacteria_Tenericutes_Mollicutes_Mollicutes_RF39_NA_NA_NA.2",
"Bacteria_Firmicutes_Negativicutes_Selenomonadales_Veillonellaceae_Anaerovibrio_NA",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_NA_NA_NA.3",
"Bacteria_Tenericutes_Mollicutes_Anaeroplasmatales_Anaeroplasmataceae_Anaeroplasma_NA",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_Lachnospiraceae_NK4A136_group_NA.6",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Ruminococcaceae_UCG.002_NA.12",
"Bacteria_Proteobacteria_Deltaproteobacteria_Desulfovibrionales_Desulfovibrionaceae_Desulfovibrio_piger",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Ruminococcus_1_NA.12",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Muribaculaceae_NA_NA.20",
"Bacteria_Cyanobacteria_Melainabacteria_Gastranaerophilales_NA_NA_NA.10",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Prevotellaceae_NK3B31_group_NA.7",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_NA_NA.108",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Ruminococcaceae_UCG.002_NA.13",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Tannerellaceae_Parabacteroides_NA.14",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Prevotella_6_NA.5",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Prevotella_9_NA.13",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Alloprevotella_NA.4",
"Bacteria_Proteobacteria_Alphaproteobacteria_Rhodospirillales_NA_NA_NA.12",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Muribaculaceae_NA_NA.21",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Ruminococcus_1_NA.14",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_Coprococcus_2_NA.4",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Bacteroidaceae_Bacteroides_NA.57",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Ruminococcus_1_NA.15",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Bacteroidaceae_Bacteroides_NA.58",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Bacteroidaceae_Bacteroides_salyersiae.1",
"Bacteria_Verrucomicrobia_Verrucomicrobiae_Verrucomicrobiales_Akkermansiaceae_Akkermansia_NA.12",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_NA_NA.133",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Alloprevotella_NA.5",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Muribaculaceae_NA_NA.24",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Muribaculaceae_NA_NA.25",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Clostridiales_vadinBB60_group_NA_NA.14",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_Anaerostipes_NA.8",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_NA_NA.154",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Tannerellaceae_Parabacteroides_NA.19",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_NA_NA.4",
"Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Prevotella_NA.7",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_NA_NA.176",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_NA_NA.87",
"Bacteria_Firmicutes_Bacilli_Lactobacillales_Lactobacillaceae_Lactobacillus_NA.13",
"Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_NA_NA.399")
rows_microbes_with_NAs <- which(rownames(data_table_4) %in% microbes_with_NAs)

}