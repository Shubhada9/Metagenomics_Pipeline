# 🧬 16S rRNA Metagenomic Analysis Pipeline

A complete R-based metagenomics pipeline covering ASV inference, diversity analysis, and taxonomic composition profiling.

## Datasets
- **DADA2 Pipeline** — Illumina MiSeq 16S rRNA (mouse gut microbiome)
- **Loue River Biofilm** — 454 pyrosequencing, 2 sites × 3 months × 3 replicates

## Experiments
| Script | Experiment | Description |
|--------|-----------|-------------|
| 01_dada2_pipeline.R | Exp 4.1 | ASV inference, chimera removal, SILVA taxonomy |
| 02_data_exploration.R | Exp 4.2 | Phyloseq object, sparsity, rarefaction curves |
| 03_alpha_diversity.R | Exp 4.3 | Richness, Chao1, Shannon, ANOVA, Kruskal-Wallis |
| 04_beta_diversity.R | Exp 4.4 | Bray-Curtis, PCoA, PERMANOVA |
| 05_taxonomic_composition.R | Exp 4.5 | Phylum pie chart + stacked bar plot |

## Key Results
- 232 ASVs generated, ~96% reads retained (DADA2)
- 18 samples, 5248 OTUs, 245,239 reads (Loue River)
- Shannon diversity significantly shaped by month (p=0.026)
- Site + month explain 69% of beta diversity variation (PERMANOVA p=0.0001)

## Requirements
```r
install.packages(c("tidyverse","vegan","patchwork","agricolae","FSA","rcompanion"))
BiocManager::install(c("phyloseq","dada2","DESeq2"))
```

## How to Run
```r
source("R/01_dada2_pipeline.R")
source("R/02_data_exploration.R")
source("R/03_alpha_diversity.R")
source("R/04_beta_diversity.R")
source("R/05_taxonomic_composition.R")
```

## Author
Shubhada Khare — MSc Bioinformatics
