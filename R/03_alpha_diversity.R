# ============================================================
# Experiment 4.3 — Alpha Diversity Analysis
# Indices: Richness, Chao1, Evenness, Shannon
# Stats: ANOVA, Tukey HSD, Kruskal-Wallis, Dunn's test
# ============================================================

library(vegan); library(phyloseq); library(tidyverse)
library(patchwork); library(agricolae); library(FSA); library(rcompanion)

data_otu <- readRDS("results/data_otu.rds")
data_grp <- readRDS("results/data_grp.rds")

# Calculate indices
data_richness <- estimateR(data_otu)
data_evenness <- diversity(data_otu) / log(specnumber(data_otu))
data_shannon  <- diversity(data_otu, index="shannon")
data_alphadiv <- cbind(data_grp, t(data_richness), data_shannon, data_evenness)
rm(data_richness)

# Tidy format
data_alphadiv_tidy <- data_alphadiv %>%
  mutate(sample_id=rownames(data_alphadiv)) %>%
  gather(key=alphadiv_index, value=obs_values, -sample_id,-site,-month,-site_month)

# Boxplots by site
P1 <- ggplot(data_alphadiv, aes(x=site,y=S.obs)) + geom_boxplot(fill=c("steelblue","tomato")) + geom_point() + labs(title="Richness",tag="A") + theme_classic()
P2 <- ggplot(data_alphadiv, aes(x=site,y=S.chao1)) + geom_boxplot(fill=c("steelblue","tomato")) + geom_point() + labs(title="Chao1",tag="B") + theme_classic()
P3 <- ggplot(data_alphadiv, aes(x=site,y=data_evenness)) + geom_boxplot(fill=c("steelblue","tomato")) + geom_point() + labs(title="Evenness",tag="C") + theme_classic()
P4 <- ggplot(data_alphadiv, aes(x=site,y=data_shannon)) + geom_boxplot(fill=c("steelblue","tomato")) + geom_point() + labs(title="Shannon",tag="D") + theme_classic()
(P1|P2)/(P3|P4)

# Richness by month x site
data_alphadiv_tidy %>%
  filter(alphadiv_index=="S.obs") %>%
  mutate(month=fct_relevel(month,"July","August","September")) %>%
  ggplot(., aes(x=month, y=obs_values)) +
  geom_boxplot(aes(fill=month)) + geom_point() +
  facet_grid(.~site) + labs(y="Richness",x="") +
  theme(axis.text.x=element_text(angle=45,hjust=1))

# ANOVA
cat("--- ANOVA: Shannon ~ site ---\n")
summary(aov(data_shannon ~ site, data=data_alphadiv))
cat("--- ANOVA: Shannon ~ month ---\n")
aov_test <- aov(data_shannon ~ month, data=data_alphadiv)
summary(aov_test)
hsd_res <- HSD.test(aov_test, "month", group=TRUE)$groups
print(hsd_res)
cat("--- Two-way ANOVA: Shannon ~ site * month ---\n")
summary(aov(data_shannon ~ site*month, data=data_alphadiv))

# Non-parametric
cat("--- Kruskal-Wallis: site ---\n")
kruskal.test(data_shannon ~ site, data=data_alphadiv)
cat("--- Kruskal-Wallis: month ---\n")
kruskal.test(data_shannon ~ month, data=data_alphadiv)
PT <- dunnTest(data_shannon ~ month, data=data_alphadiv, method="bh")
cldList(comparison=PT$res$Comparison, p.value=PT$res$P.adj, threshold=0.05)

write.csv(data_alphadiv, "results/alpha_diversity_indices.csv")
