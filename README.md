# PEER_factors
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

# Citation

Angli Xue, Seyhan Yazar, Drew Neavin, Joseph E. Powell. Pitfalls and opportunities for applying latent variables in single-cell eQTL analyses. _provisionally accepted_. Genome Biology. 2023. [BioRxiv](https://www.biorxiv.org/content/10.1101/2022.08.02.502566v1)

For questions, please email us at Angli Xue (a.xue@garvan.org.au) or Joseph E. Powell (j.powell@garvan.org.au)
