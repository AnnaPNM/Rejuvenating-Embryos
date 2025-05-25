# Rejuvenating-Embryos

## Embryogenesis

For embryogenesis we needed to have data about as many time points as possible. We selected a reference study https://www.nature.com/articles/s41592-024-02511-3 where authors have already collected a lot of RNA-seq embryonic datasets, then merged them and corrected the butch-effect. As a result, they created a scanpy-readable h5ad object. This source contains data about both mouse and human. In this way, these single-cell datasets (both mouse and human) were selected for embryogenesis trajectories building.


## Reprogramming

A list of reprogramming datasets was compiled, four of which were selected for further analysis. For each dataset, raw count matrices were obtained if not publicly accessible, metadata were assembled, a variance‐stabilizing transformation (VST) was applied, and Spearman’s correlation coefficients were computed for each gene and hallmark. The code implementing these transformations is provided in the Trajectories/Reprogramming folder.


## Aging

To integrate multiple bulk RNA-seq datasets related to aging in Mus musculus, we merged raw count data and applied ComBat-seq, from sva R-package to correct for batch effects. We then used DESeq2 R-package with a design formula including Day and Tissue to account for variation due to those factors, using these counts later in a likelihood ratio test (LRT) to identify time-dependent genes. We also applied variance stabilizing transformation (VST) from DESeq2 to the ComBat-seq-corrected counts, and used the resulting data throughout the rest of the analysis - including expression trajectory plotting and co-expression module detection. For aging in Homo Sapiens, we needed data for the different ages. We chose the Tabula Sapiens dataset (https://figshare.com/articles/dataset/Tabula_Sapiens_v2/27921984) because it contains early ages of different cell types from 11 male and 13 female donors. They  balance and assign cell types from each tissue compartment and optimally mix high-quality plate-seq data and high-volume droplet-based data to provide a broad and deep benchmark atlas. As a result, they created a scanpy-readable object called h5ad.
