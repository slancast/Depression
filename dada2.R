#!/user/bin/R
#from https://benjjneb.github.io/dada2/tutorial.html

library(dada2); packageVersion("dada2")
path <- "/Users/slancast/Box/School\ for\ The\ Work\ Data\ --\ Sam/Microbiome/All_fastqs" 
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnFs <- fnFs[lapply(fnFs,function(x) length(grep(".gz",x,value=FALSE))) == 0] #removing the .gz files
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
fnRs <- fnRs[lapply(fnRs,function(x) length(grep(".gz",x,value=FALSE))) == 0] #removing the .gz files

sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1) #Creating a list of names cutting off the last characters from the file name

plotQualityProfile(fnFs[7:8])

plotQualityProfile(fnRs[3:4])

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#```{r assign filenames for de-primered files}
# Place files without primers in clean/ subdirectory
noprimerFs <- file.path(path, "clean", paste0(sample.names, "_F_noprimer.fastq.gz"))
noprimerRs <- file.path(path, "clean", paste0(sample.names, "_R_noprimer.fastq.gz"))
names(noprimerFs) <- sample.names
names(noprimerFs) <- sample.names
#```

#```{r remove primers}
# forward primer: 515F: GTGCCAGCMGCCGCGGTAA
# reverse primer: 806R: GGACTACHVGGGTWTCTAAT
FWD_PRIMER <- "GTGCCAGCMGCCGCGGTAA"
REV_PRIMER <- "GGACTACHVGGGTWTCTAAT"
FWD_PRIMER_LEN = nchar(FWD_PRIMER)
REV_PRIMER_LEN = nchar(REV_PRIMER)
trimLeft = c(FWD_PRIMER_LEN, REV_PRIMER_LEN)

out <- filterAndTrim(fwd = fnFs[1:20], filt = noprimerFs[1:20], rev = NULL, filt.rev = NULL,
                     trimLeft = FWD_PRIMER_LEN, maxN=0, maxEE=2, truncLen=150,
                     rm.phix=TRUE,compress=TRUE, multithread=TRUE)

# On Windows set multithread=FALSE
head(out)
#```

#```{r inspect read quality}
# plot quality of the forward reads
plotQualityProfile(noprimerFs[1:3])

# plot quality of the reverse reads
plotQualityProfile(noprimerRs[1:3])
#```


#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(140,140), #the default truncLen didn't work for this project
                     #maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, 
                     #compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

#head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE) #Ariel said Holms said to only do forwared reads. If so we take this out.

#seqtab <- makeSequenceTable(mergers) #Holy cow the merged paris don't work at all. Not clear to me exactly why. Looks like something happens when merging.
seqtab <- makeSequenceTable(dadaFs)#This is what I think will need to do with only forward reads

dim(seqtab)

table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#The axonomy() function only works wtih sequences that are above 50bp, so I'll implement the following code
reads_over50 <- which(nchar(colnames(seqtab.nochim)) >= 50)
seqtab.nochim <- seqtab.nochim[,reads_over50]

taxa <- assignTaxonomy(seqtab.nochim, paste0(path,"/silva_nr_v132_train_set.fa.gz"), multithread=TRUE) #, tryRC=TRUE

taxa <- addSpecies(taxa, paste0(path,"/silva_species_assignment_v132.fa.gz"))
taxa_t <- t(taxa)
count_output <- rbind(taxa_t, seqtab.nochim)
write.table(count_output, file = paste0("/Users/slancast/Desktop/depression_dada2.txt"),sep="\t", quote=F, col.names = F)


taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

unqs.mock <- seqtab.nochim["224325473",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

library(ggplot2); packageVersion("ggplot2")
library(Biostrings); packageVersion("Biostrings")
library(phyloseq); packageVersion("phyloseq")

theme_set(theme_bw())

# This is them adding metadta.
# 
# samples.out <- rownames(seqtab.nochim)
# subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
# gender <- substr(subject,1,1)
# subject <- substr(subject,2,999)
# day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
# samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
# samdf$When <- "Early"
# samdf$When[samdf$Day>100] <- "Late"
# rownames(samdf) <- samples.out


ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

plot_richness(ps, measures=c("Shannon"))



ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, title="Bray NMDS")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20,  fill="Family") # + facet_wrap(~When, scales="free_x") x="Day",

