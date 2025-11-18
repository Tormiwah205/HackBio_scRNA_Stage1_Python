# Literature Review: Single-Cell RNA-seq Analysis

This is a review of the assigned essay on "Best Practices in Single-Cell RNA-seq Analysis: A Tutorial".

## Best Practices in Single-Cell RNA-seq Analysis: A Tutorial

Single-cell RNA sequencing (scRNA-seq) has fundamentally advanced the study of biological systems, enabling the analysis of gene expression at an unprecedented resolution and leading to the discovery of previously obscured cellular populations. However, the relative immaturity of the field and this proliferation of methods present significant challenges to establishing standardized scRNA-seq analysis workflows. To address this, current best practices outline a comprehensive workflow that progresses from raw count matrices through downstream analysis. 

The initial stage, pre-processing, involves several steps: quality control (QC), normalization, data correction, feature selection, and dimensionality reduction. Cell QC is crucial, requiring the joint examination of key covariates, to identify and filter out outliers like dying cells or doublets. For normalization, which aims to correct for variability in captured mRNA molecules due to sampling effects, the scran pooling-based size factor estimation method is highly recommended for non-full-length datasets, consistently demonstrating top performance in comparison studies. Subsequently, data correction addresses unwanted variability; for straightforward batch effects within a single experiment, the linear method ComBat is recommended as a strong performer. A critical recommendation is to perform statistical testing on "measured data" rather than on "corrected data," which is better reserved for visual comparison. For efficiency in downstream analysis, dimensionality reduction is performed. For visualization, non-linear methods are preferred; UMAP (Uniform Approximation and Projection) is regarded as best practice for exploratory data visualization due to its speed and ability to scale to large datasets.

Downstream analysis focuses on biological inference through cell-level and gene-level approaches. Cell clustering, often achieved using Louvain community detection on a K-Nearest Neighbour (KNN) graph, has been shown to be robust for scRNA-seq data. Trajectory inference (TI) models continuous biological processes; the Slingshot method is highly effective for reconstructing simpler linear or bifurcating trajectories. Finally, differential expression (DE) testing is performed, typically within cell-identity clusters. Recommended tools include the single-cell specific hurdle model MAST or bulk methods adapted for single-cell data.

### Reference

Luecken, M. D., & Theis, F. J. (2019). Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology, 15*(e8746). https://doi.org/10.15252/msb.20188746
