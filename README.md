# PEER_factors
This repository contains the analysis code pipeline to generate PEER factors from pseudo-bulk data and perform eQTL association analysis as part of the manuscript "**Pitfalls and opportunities for applying PEER factors in single-cell eQTL analyses**"

Scripts are listed by the order in the methods section of the manuscript:

1. Extract the whole OneK1K dataset from .RDS and subgroup into 14 cell types
2. Generate the pseudo-bulk mean matrix
3. Generate PEER factors with 13 QC options

	a. Extra information of runtime and nr of iterations
	
	b. Make new covariate files
4. Run sensitivity test by MatrixeQTL
	
	a. Merge results
	
	b. Summarize and nr of eQTL and eGenes
5. Down-sampling analysis



All code is also available on Github: https://github.com/anglixue/PEER_factors/
