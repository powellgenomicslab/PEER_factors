# Pitfalls and opportunities for applying latent variables in single-cell eQTL analyses
This repository contains the analysis code pipeline to generate PEER factors from pseudo-bulk data and perform eQTL association analysis as part of the manuscript "**Pitfalls and opportunities for applying latent variables in single-cell eQTL analyses**"

Scripts are listed by the order in the methods section of the manuscript:

1. Extract the whole OneK1K dataset from .RDS and subgroup into 14 cell types
2. Generate the pseudo-bulk mean matrix
3. Generate PEER factors (PFs) with 13 QC options
    + Extra information of runtime and nr of iterations
    + Make new covariate files
4. Run sensitivity test by MatrixeQTL
    + Merge results
    + Summarize and nr of eQTL and eGenes
5. Down-sampling analysis
6. Principal component analysis (PCA)
    + Generate PCs by [PCAForQTL](https://github.com/heatherjzhou/PCAForQTL)
    + Run eQTL sensitivity test adjusting PCs 0-50 (similar to PFs)
7. Main figures
8. Supplementary figures and tables

Folder v1.0 contains the older version of the scripts.

All code is also available on Angli's personal Github: https://github.com/anglixue/PEER_factors/ and Zenodo: https://doi.org/10.5281/zenodo.7513270

# OneK1K Seurat object

We have submitted the OneK1K Seurat object that was used in Yazar et al. (Science. 2022) to the cellxgene. However, the platform only accepts GRCh38 but in Yazar et al we used GRCh37. So we have to realign the data to fulfill the website requirement, which explains the slight difference in the number of cells between the published version and cellxgene version.

Also, we suggest removing two individuals (“88_88” and “966_967”) due to the extremely low number of cells and abnormal cell composition. The exact steps can be found on line #18-22 on this [Github page](https://github.com/powellgenomicslab/PEER_factors/blob/main/1-Extract_datasets/Extract_RDS_all_cell_types.R).

If you need the updated Seurat object of OneK1K, please email Angli Xue (a.xue@garvan.org.au).

# Citation

Angli Xue, Seyhan Yazar, Drew Neavin, Joseph E. Powell. Pitfalls and opportunities for applying latent variables in single-cell eQTL analyses. _Genome Biology_. 2023. [[Full text](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02873-5)]

For questions, please email us at Angli Xue (a.xue@garvan.org.au) or Joseph E. Powell (j.powell@garvan.org.au).
