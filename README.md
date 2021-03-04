# TCGA data mining for extablishing the relationship between sub-branches of Wnt signalling in colorectal cancer (CRC)

Project completed during my BSc in Biochemistry at the University of Southampton, 2018-2019.

## Data extraction
* produce a library of Wnt signalling components and their associated genes of proteins
* search for CRC gene expression data (e.g. RNAseq) in publicly-available datasets (e.g. TCGA)

## Data pre-processing
* ensure a fair data comparison (data was normalised using RSEM - RPKM modelled to account for isoform abundances)
* address the problem of missing data

## Data analysis
* uses ***gplots**, ***RColorBrewer*** and ***dunn.test*** R packages
1. **EXPLORATORY DATA ANALYSIS**
   * check distribution of data
3. **UNSUPERVISED MACHINE LEARNING**
   * hierarchical clustering (no prior information was known about the groups)
   * Pearson correlation heatmap (for computing distances using complete linkage)
4. **CLUSTER ANALYSIS FOR SIGNIFICANT GROUPS**
   * Kruskal-Wallis test for comparing groups
   * post-hoc Dunn's test with Benjamini-Hochberg p-value correction for finding which group was different
   * hypergeometric test for the probabilities of randomly picking the samples from a certain category
5. **FURTHER ANALYSIS** 
   * gene set enrichment (PANTHER)
   * PPIs network (STRING)
