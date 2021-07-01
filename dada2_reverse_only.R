seqtab <- makeSequenceTable(dadaRs)#This is what I think will need to do with only reverse reads

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

#The assignTaxonomy() function only works wtih sequences that are above 50bp, so I'll implement the following code
reads_over50 <- which(nchar(colnames(seqtab.nochim)) >= 50)

taxa <- assignTaxonomy(colnames(seqtab.nochim)[reads_over50], paste0(path,"/silva_nr_v132_train_set.fa.gz"), multithread=TRUE, tryRC=TRUE)

taxa <- addSpecies(taxa, paste0(path,"/silva_species_assignment_v132.fa.gz"))

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
               sample_data(samdf), 
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

