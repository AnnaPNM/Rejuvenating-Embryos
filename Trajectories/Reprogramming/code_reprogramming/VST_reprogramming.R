if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "pheatmap", "ggplot2"), ask = FALSE)

library(DESeq2)
library(pheatmap)
library(ggplot2)

count_data <- read.csv("../data_example_reprogramming/example_raw_counts.csv", row.names = 1)
# data_example: raw count matrix where row indices are gene names and column names are samples

count_matrix <- t(as.matrix(count_data))

metadata <- read.csv("../data_example_reprogramming/example_metadata.csv",
                     row.names = 1)
# metadata is manually curated and should include columns 'days_of_reprogramming' and 'cell_type'. 
# in 'days_of_reprogramming', each sample should have numeric values, except for iPSC samples, which retain the string 'iPSC'.

metadata$days_of_reprogramming <- as.factor(metadata$days_of_reprogramming)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ days_of_reprogramming)

keep_genes <- rowSums(counts(dds)) >= 10
dds <- dds[keep_genes, ]

dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = TRUE)

vsd <- vst(dds, blind = TRUE)
vst_matrix <- assay(vsd)

pcaData <- plotPCA(vsd, intgroup = "days_of_reprogramming", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = days_of_reprogramming)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% of variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% of variance")) +
  ggtitle("PCA after VST")
ggsave("../output_reprogramming/PCA_after_VST_reprogramming.png", width = 6, height = 5, dpi = 300)

topVarGenes <- head(order(rowVars(vst_matrix), decreasing = TRUE), 50)
png("../output_reprogramming/heatmap_top50_genes_example.png", width = 1000, height = 800, res = 150)
pheatmap(vst_matrix[topVarGenes, ],
         annotation_col = metadata,
         scale = "row",
         clustering_distance_rows = "correlation",
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()

write.csv(assay(vsd), gzfile("../output_reprogramming/example_vst.csv.gz"))


