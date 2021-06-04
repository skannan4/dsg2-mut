# Electrophysiological defects in a novel patient-derived stem cell model of desmoglein-2 mutant ARVC

### Introduction
Arrhythmogenic right ventricular cardiomyopathy (ARVC) is a progressive, inheritable heart condition that leads to fibro-fatty scarring of the myocardium and causes ventricular arrhythmias and sudden cardiac death. There has been significant interest in understanding the pathology of ARVC, but given the generally poor animal models, and the fact that the disease process likely begins before the clinical presentation, there has been an interest in developing better models of human disease. Pluripotent stem-cell derived cardiomyocytes (PSC-CMs) may serve as one such model. In this study, Rob Hawthorne et al. (Les Tung lab) generated a patient-derived DSG2-mutant line to study ARVC pathogenesis. As part of this study, we performed bulk RNA-seq to look at transcriptomic differences between DSG2-mutant iPSC-CMs compared to iPSC-CMs generated from an age- and sex-matched control.

### Data
All of the relevant files, including the counts table, phenotype table, and a workspace with preloaded objects, can be found on [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165300) (this isn't the correct GEO, correct when uploaded). The one additional file that can be found here is goi.txt, which just contains genes of interest for plotting (in this particular case, they were genes drawn from the NFKB pathway).

### Dependencies
Most of the libraries used in our codebase can be found from CRAN or Bioconductor. However, we additionally make use of the SingleCellNet package from the Cahan lab. Please see [their github](https://github.com/pcahan1/singleCellNet) for instructions on how to install SingleCellNet.
