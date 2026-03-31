# ============================================================
# Beta Diversity Analysis
# Bray-Curtis dissimilarity, PCoA, PERMANOVA
# ============================================================

library(vegan); library(phyloseq); library(tidyverse)

data_phylo <- readRDS("results/data_phylo.rds")
data_grp   <- readRDS("results/data_grp.rds")
TAX <- tax_table(data_phylo)
SAM <- sample_data(data_phylo)

# Filter low-abundance OTUs
data_phylo_filt <- filter_taxa(data_phylo,
  function(x) sum(x>2) > (0.11*length(x)), TRUE)

# Rarefy
set.seed(1782)
OTU_filt_rar <- rarefy_even_depth(otu_table(data_phylo_filt), rngseed=TRUE, replace=FALSE)
data_otu_filt_rar  <- data.frame(otu_table(OTU_filt_rar))
data_phylo_filt_rar <- phyloseq(OTU_filt_rar, TAX, SAM)

# Bray-Curtis matrix
dist_bc <- as.matrix(vegdist(data_otu_filt_rar, method="bray"))
print(round(dist_bc[1:5,1:5], 3))

# PCoA by site
pcoa_bc <- ordinate(data_phylo_filt_rar, "PCoA", "bray")
plot_ordination(data_phylo_filt_rar, pcoa_bc, color="site") +
  geom_point(size=4) + theme_bw() +
  labs(title="PCoA — Bray-Curtis (by site)")

# PCoA by site + month
plot_ordination(data_phylo_filt_rar, pcoa_bc, color="site", shape="month") +
  geom_point(size=4) + theme_bw() +
  labs(title="PCoA — Bray-Curtis (site + month)")

# PERMANOVA
cat("--- PERMANOVA: site ---\n")
adonis2(data_otu_filt_rar~site, data=data_grp, permutations=9999, method="bray")
cat("--- PERMANOVA: site * month ---\n")
adonis2(data_otu_filt_rar~site*month, data=data_grp, permutations=9999, method="bray")

saveRDS(data_otu_filt_rar,   "results/data_otu_filt_rar.rds")
saveRDS(data_phylo_filt_rar, "results/data_phylo_filt_rar.rds")
