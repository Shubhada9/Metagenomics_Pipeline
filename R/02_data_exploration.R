# ============================================================
# Experiment 4.2 — Data Exploration & Quality Assessment
# Loue River Epilithic Biofilm 16S rRNA data
# ============================================================

library(vegan); library(phyloseq); library(tidyverse)

# Import data — UPDATE PATHS
data_otu  <- read.table("data/data_loue_16S_nonnorm.txt", header=TRUE)
data_grp  <- read.table("data/data_loue_16S_nonnorm_grp.txt", header=TRUE, stringsAsFactors=TRUE)
data_taxo <- read.table("data/data_loue_16S_nonnorm_taxo.txt", header=TRUE)

# Build phyloseq object
OTU <- otu_table(as.matrix(data_otu), taxa_are_rows=FALSE)
SAM <- sample_data(data_grp)
TAX <- tax_table(as.matrix(data_taxo))
data_phylo <- phyloseq(OTU, TAX, SAM)
data_phylo

# Basic stats
nb_samples <- dim(data_otu)[1]
nb_var     <- dim(data_otu)[2]
cat("Samples:", nb_samples, "| OTUs:", nb_var, "| Total reads:", sum(data_otu), "\n")

# Sparsity
cat("Zeros:", sum(data_otu==0), "| Sparsity:", round(sum(data_otu==0)/(nb_var*nb_samples)*100,2), "%\n")

# Plots
hist(as.matrix(data_otu), max(data_otu), right=FALSE, las=1,
     xlab="Occurrence value", ylab="Frequency", main="OTU Count Frequency")

non_zero <- sapply(1:nb_var, function(i) sum(data_otu[,i] != 0))
plot(non_zero, xlab="OTU", ylab="Non-zero count", main="OTU Prevalence", pch=20)

rarecurve(data_otu, step=100, cex=0.75, las=1)

sum_seq <- rowSums(data_otu)
plot(sum_seq, ylim=c(0,25000), main="Sequencing Depth per Sample", pch=19)
cat("Min:", min(sum_seq), "| Max:", max(sum_seq), "\n")

# Save
saveRDS(data_phylo, "results/data_phylo.rds")
saveRDS(data_otu,   "results/data_otu.rds")
saveRDS(data_grp,   "results/data_grp.rds")
saveRDS(data_taxo,  "results/data_taxo.rds")
