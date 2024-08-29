# aptpro
 This package provides tools for quality control of single-cell nucleic acid aptamers, gRNA, and mRNA multi-omics sequencing. It also includes functions for defining cell gRNA identities and predicting potential proteins that may bind to nucleic acid aptamers families.
## Installation
You can install this R package using the following methods:

Use the `devtools` package to install the latest version directly from GitHub:

'''
#install.packages("devtools")
devtools::install_github("fuyboo/aptpro")

## Example
### Step1：Nucleic Acid Aptamer Clustering and Family Classification
  Aptamer binding interactions are primarily driven by conserved regions, with minimal differences observed in single aptamer interactions with various proteins. To address this, we use the clustering strategy proposed in Smart-Aptamer, which focuses on selecting the most abundant sequences for family clustering.

'''
python ./code/smart_cluster.py  -t 35 -i 0.7 -e 0.05 -o ./lgy/data_3/motif/test1w
