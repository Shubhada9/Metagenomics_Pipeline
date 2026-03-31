# ============================================================
# Experiment 4.1 — DADA2 Pipeline
# High-resolution ASV inference from Illumina MiSeq 16S rRNA
# ============================================================

library(dada2)

# Set path — UPDATE THIS
path <- "path/to/MiSeq_SOP"
list.files(path)

# Forward and reverse fastq files
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names=TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names=TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Filter and trim
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
        truncLen=c(240,160), maxN=0, maxEE=c(2,2),
        truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE)
head(out)

# Learn error rates
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)
plotErrors(errF, nominalQ=TRUE)

# Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Sequence table
seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                  multithread=FALSE, verbose=TRUE)
cat("Fraction retained:", sum(seqtab.nochim)/sum(seqtab), "\n")

# Track reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs,getN), sapply(dadaRs,getN),
               sapply(mergers,getN), rowSums(seqtab.nochim))
colnames(track) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
rownames(track) <- sample.names
head(track)

# Assign taxonomy — UPDATE PATH
taxa <- assignTaxonomy(seqtab.nochim,
         "path/to/silva_nr99_v138.1_train_set.fa.gz", multithread=FALSE)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(seqtab.nochim, "results/seqtab_nochim.rds")
saveRDS(taxa, "results/taxa.rds")
