# Genetic effects of tissue-specific enhancers in schizophrenia and hypertrophic cardiomyopathy
## Doctoral thesis by Emanuele Felice Osimo for Imperial College London
The link to the thesis will be included once publicly deposited.
Contact me at eosimo at ic.ac.uk

### Chapter Chapter 2 - Schizophrenia and HCM heritability enrichment in tissue-specific enhancers

This repository contains the code for work contained in Chapter 2. This is divided into three files, which correspond to work packages:

- _1_enhancer_list_generation.Rmd_ is R code for generating tissue-expressed lists of enhancers based on AR+C and FANTOM5 annotations.
- _2_LDSC.sh_ is bash code for running LDSC analyses on these lists.
- _3_LDSC_plot_results.R_ is R code for plotting LDSC results for publication.
- The _E_P_files_ folder contains genome-wide annotations for all FANTOM5 enhancers, with coordinates extended by 100bp at each end. Neural- and cardiac-expressed enhancers, with annotations from the AR+C, will be provided upon request to the Author.