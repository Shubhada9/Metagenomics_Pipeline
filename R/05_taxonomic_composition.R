# ============================================================
# Taxonomic Community Composition
# Global pie chart + per-treatment stacked bar plot
# ============================================================

library(vegan); library(phyloseq); library(tidyverse)

data_phylo        <- readRDS("results/data_phylo.rds")
data_grp          <- readRDS("results/data_grp.rds")
data_taxo         <- readRDS("results/data_taxo.rds")
data_otu_filt_rar <- readRDS("results/data_otu_filt_rar.rds")
TAX <- tax_table(data_phylo)
SAM <- sample_data(data_phylo)

# Filter taxonomy to match rarefied OTUs
data_taxo_filt_rar <- data_taxo[rownames(data_taxo) %in% colnames(data_otu_filt_rar),]
data_taxo_filt_rar$OTU_id <- rownames(data_taxo_filt_rar)

# Merge OTU + metadata
data_grp_temp <- data_grp; data_grp_temp$sample_id <- rownames(data_grp_temp)
data_otu_temp <- data_otu_filt_rar; data_otu_temp$sample_id <- rownames(data_otu_temp)
data_otu_grp  <- inner_join(data_grp_temp, data_otu_temp, by="sample_id")
rm(data_grp_temp, data_otu_temp)

# Aggregate by treatment
data_agg <- aggregate(data_otu_grp[,5:ncol(data_otu_grp)],
              by=list(site_month=data_otu_grp$site_month), sum)

# Transpose and add totals
data_agg_t <- as.data.frame(t(as.matrix(data_agg[,2:ncol(data_agg)])))
colnames(data_agg_t) <- data_agg[,1]
data_agg_t$total_counts <- rowSums(data_agg_t)
data_agg_t$OTU_id <- rownames(data_agg_t)

data_otu_t <- as.data.frame(t(as.matrix(data_otu_filt_rar)))
data_otu_t$OTU_id <- rownames(data_otu_t)
data_merged <- inner_join(data_agg_t, data_otu_t, by="OTU_id")
data_final  <- inner_join(data_taxo_filt_rar, data_merged, by="OTU_id")

# Global pie chart
abundant_com <- data_final %>%
  group_by(LCA_simplified) %>%
  summarize(total=sum(total_counts)) %>%
  mutate(pct=total/sum(total)*100) %>%
  filter(pct >= 1) %>% select(LCA_simplified, pct)
rare_com <- tibble(LCA_simplified="Rare_Phyla", pct=100-sum(abundant_com$pct))
global_com <- bind_rows(abundant_com, rare_com)

ggplot(global_com, aes(x="", y=pct, fill=LCA_simplified)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + theme_void() +
  labs(title="Global Bacterial Community Composition") +
  theme(legend.position="bottom", legend.title=element_blank())

# Stacked bar plot per treatment
treatment_cols <- c("Cleron_July","Cleron_August","Cleron_September",
                    "Parcey_July","Parcey_August","Parcey_September")
com_phylum <- aggregate(data_final[,treatment_cols],
               by=list(Phylum=data_final$Phylum), sum) %>%
               gather(key=treatment, value=counts, -Phylum)

com_phylum %>%
  mutate(treatment=fct_relevel(treatment,
    "Cleron_July","Cleron_August","Cleron_September",
    "Parcey_July","Parcey_August","Parcey_September")) %>%
  ggplot(., aes(x=treatment, y=counts, fill=Phylum)) +
  geom_bar(position="fill", stat="identity", color="white") +
  scale_y_continuous(labels=scales::percent_format()) +
  labs(title="Bacterial Composition at Phylum Level", x="", y="Proportion") +
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1),
                     legend.title=element_blank())

write.csv(global_com, "results/global_composition.csv")
