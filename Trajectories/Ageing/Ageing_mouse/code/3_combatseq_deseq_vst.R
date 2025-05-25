library(sva)
library(readr)

library(BiocParallel)
register(MulticoreParam(2))

counts <- read_csv("counts_matrix2.csv")
meta <- read_csv("metadata.csv")

counts <- as.data.frame(counts)
rownames(counts) <- counts[[1]]
counts <- counts[,-1]

counts <- counts[rowSums (counts) > 0,]

combat_counts <- ComBat_seq(counts = as.matrix(counts),
                            batch = meta$Dataset)

write.csv(combat_counts,'combatseq_corrected_counts.csv', row.names = TRUE)

#DESeq2

combat_counts <- read_csv("combatseq_corrected_counts.csv")

library(DESeq2)

meta$Tissue <- factor(meta$Tissue)
meta$Day <- factor(meta$Day)

dds <- DESeqDataSetFromMatrix(countData = combat_counts, colData = meta, design = ~ Day + Tissue)

dds <- estimateSizeFactors(dds, type = "poscounts")

dds <- DESeq(dds, parallel = TRUE)

saveRDS(dds, file = "dds_preLRT.rds")

norm_counts <- counts(dds, normalized = TRUE) #save deseq normalised

write.csv(norm_counts,'DESeq_combat_corrected_counts.csv')

#Test for global differences across all timepoints
# (If Day has many levels (e.g., day2, day10, ..., day900), 
# and you want genes that change across time (not just between two specific days), 
# use the likelihood ratio test (LRT):)

dds <- readRDS("dds_preLRT.rds")

dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ Tissue, parallel = TRUE)
res <- results(dds_lrt)

res <- res[order(res$padj), ]
res <- res[!is.na(res$padj), ]

res_sig <- res[res$padj < 0.05, ]
# Sort and extract DE genes
top_up <- head(res_sig[order(-res_sig$log2FoldChange), ], 300)
top_down <- head(res_sig[order(res_sig$log2FoldChange), ], 300)

top_genes <- c(rownames(top_up), rownames(top_down))

# Save
write.csv(as.data.frame(res), file = "deseq2_results.csv")
writeLines(top_genes, "Ageing_top_genes.txt")

#VST 

combat_counts <- read_csv("combatseq_corrected_counts.csv", row.names = 1)
meta <- read_csv("metadata.csv")


combat_counts[combat_counts == 0] <- 1
dds <- DESeqDataSetFromMatrix(countData = combat_counts, colData = meta, design = ~ Day + Tissue)

vsd <- vst(dds, blind=TRUE)
vst_mat <- assay(vsd)
write.csv(as.data.frame(vst_mat), file = "vst_normalized_counts.csv", row.names = TRUE)

plotPCA(vsd, intgroup = "Tissue")

