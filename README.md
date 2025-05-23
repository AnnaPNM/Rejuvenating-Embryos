# Rejuvenating-Embryos

## Reprogramming

Raw count matrices were generated using the nf-core pipelines rnaseq and scrnaseq with default parameters. In the final count matrices, rows are gene identifiers (ENSEMBL IDs) and columns are sample names. Metadata are manually curated and must include the following columns: days_of_reprogramming – numeric values for each sample, except iPSC samples which retain the string "iPSC"; and cell_type.

1. VST_reprogramming.R

Performs DESeq2-based normalization (variance stabilizing transform), PCA, and heatmap on raw count data.

2. trends_aggregation_reprogramming.py

Aggregates expression trajectories for a given subprocess via Z-score normalization and PCA, then plots the first principal component over time.

```
`data`: DataFrame of expression values (samples × genes)
`autophagy`: DataFrame of gene annotations for the "autophagy" subprocess
`metadata`: DataFrame with 'days_of_reprogramming'
plot_subprocess_trends(
    data=data,
    genes=autophagy,
    metadata=metadata,
    ensembl_column="ENSEMBL",
    name_subpocess="Autophagy"
)
```

3. spearman_reprogramming.py

compute_hallmark_spearman

```
df_genes = compute_hallmark_spearman(
    metadata=metadata,
    expr_data=data,
    all_hallmarks=all_hallmarks
)
```

What it does:

Filters Hallmark genes present in expr_data.

Computes Spearman’s ρ and raw p-value vs. days_of_reprogramming.

Applies FDR correction and annotates with gene symbols.

Returns:
DataFrame with columns ['Hallmark', 'Subprocess', 'GeneSymbol', 'SpearmanR', 'raw_pval'], sorted by Hallmark and Subprocess.

compute_subprocess_spearman

```
df_subproc = compute_subprocess_spearman(
    metadata=metadata,
    expr_data=data,
    all_hallmarks=all_hallmarks,
    fdr_alpha=0.05,
    min_genes=3
)
```

What it does:

Computes gene-level Spearman’s ρ and p-values as above.

Performs FDR correction and flags significant genes.

Aggregates by (Subprocess, Hallmark) pairs to compute: median SpearmanR, median raw p-value, median FDR, total number of genes (N_Genes), and number of significant genes (N_Significant).

Returns:
DataFrame with columns ['Subprocess', 'Hallmark', 'MedianSpearmanR', 'MedianPValue', 'MedianFDR', 'N_Genes', 'N_Significant'], sorted by Hallmark and Subprocess.
