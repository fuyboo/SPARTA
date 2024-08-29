# aptpro
 This package provides tools for quality control of single-cell nucleic acid aptamers, gRNA, and mRNA multi-omics sequencing. It also includes functions for defining cell gRNA identities and predicting potential proteins that may bind to nucleic acid aptamers families.
## Installation
You can install this R package using the following methods:

Use the `devtools` package to install the latest version directly from GitHub:

```
#install.packages("devtools")
devtools::install_github("fuyboo/aptpro")
```

## Example
### Step1：Nucleic Acid Aptamer Clustering and Family Classification
  Aptamer binding interactions are primarily driven by conserved regions, with minimal differences observed in single aptamer interactions with various proteins. To address this, we use the clustering strategy proposed in Smart-Aptamer, which focuses on selecting the most abundant sequences for family clustering.

```
python ./code/smart_cluster.py  -t 35 -i 0.7 -e 0.05 -o ./lgy/data_3/motif/test1w
```

### Step2: Cell Quality Control
  Before performing cell quality control, ensure that you have a Seurat object that includes three essential components: mRNA, aptamer, and motif.

```
SUM159<-cell_quality (SUM159,
                      count_threshold = 100,
                      feature_threshold = 200,
                      percent_mt_threshold = 10,
                      assay = "sgRNA",
                      save_path = NULL)
```

### Step3: Define Cell gRNA Identity
  In this step, you will define the cell gRNA identity based on the quality and enrichment ratios of the cell gRNA, which were assessed in Step 1. This process involves setting thresholds to categorize cell gRNA effectively.

```
SUM159<-cell_gRNA_identity(SUM159,
                           min_count = 200,
                           min_ratio = 0.7)
```

### Step4：Predict Aptamer Family Protein Binding
  In this step, you will predict which proteins are likely bound by the aptamer families. This involves calculating the differential matrix based on the median difference between target cells and NC (negative control) cells, filtering out low-difference families and confusing targets, and using a Gaussian Mixture Model (GMM) to refine the predictions.

```
predict_result<-predict_apt_pro(SUM159,
                                assay = "motif1w",
                                top_n = 20,
                                save_path = NULL)
```

### Step5: Visualize Prediction Results
  Visualization is crucial for interpreting and presenting the results of your aptamer-protein binding predictions. In this step, you'll generate plots to help understand the data and findings from the previous analysis.

```
visualize_aptamer_difference(predict_result,'Clust-1')

```


