# multi-omics_integration
R scripts for multi-omics statistical integration


## dependencies 
R packages :
- Mixomics v3.9 http://www.bioconductor.org/packages/release/bioc/html/mixOmics.html
- Vegan v2.5-5 https://cran.r-project.org/web/packages/vegan/index.html

## overview OTU-metabolites_toydatasets
Script for the statistical integration of OTU taxonomic counts and metabolomic dataset performed on the same samples. 
After a first exploration using PCA and procrustes analysis, a series of PLS and PLS-DA are performed on the dataset. Finally a multilevel PLS and PLS-DA was performed on the datasets to take into account the time-series.
