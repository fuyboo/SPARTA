# SPARTA
This package offers tools for SPARK-seq data analysis. It includes three main steps: 1) quality control and cell classification, 2) aptamer-target interaction prediction, and 3) data visualization.

Additionally, SPARTA includes two deep learning modules: one for aptamer binding potential prediction and another for de novo aptamer generation.
 




![](picture/model_diagram%20.png "annotation")




## Installation
First, install the necessary dependencies for aptamer family clustering, which are MCL (https://micans.org/mcl/) and BLAST (https://blast.ncbi.nlm.nih.gov/).

```
conda install bioconda::mcl
conda install biopython
```

You can install this R package using the following method:
Use the `devtools` package to install the latest version directly from GitHub:

```
#install.packages("devtools")
devtools::install_github("fuyboo/aptpro")
```

## Input data preparation


### Raw data preparation
The process from raw data to the generation of mRNA, aptamer, and sgRNA expression matrices can be referred to in `./raw_process/raw_process.pdf`.


### Aptamer Family Classification
Based on the aptamer results from the previous step, a certain number of aptamers can be selected for family analysis. For example, we can select the top 10,000 most abundant sequences for family clustering, as referenced in `./data/input/uniq_aptamer.fasta`.


Classify aptamer sequences based on their similarities using the BLAST-vs-BLAST and MCL strategy.
-t threads -i inflation value for mcl algorithm  -e pvalue_threshold -o output_directory

```

python ./aptamer_family_analysis/smart_cluster.py  -t 35 -i 0.7 -e 0.05 -o ./lgy/data_3/motif/test1w

```

Finally, the corresponding family groups of the aptamers are saved in a file such as ./data/output/Aptamer_family.csv.


| name  | seq | seq | 
| ------------- | ------------- | ------------- |
| Apt-1  | TTTCGGCGGGTGAATATCCAACTGGTCCGTCCCTTGGGATCTTTGT  | Clust-5  |
| Apt-2  | GGTTTGCTGAGGTGGGCGTCGTTGAATGTTAGTTCGGGAATACTTG  | Clust-3  |
| Apt-3  | GGCTCCTCTTAGGGGCTGTGACCGGCGGGCGGGAATGTAGCAGGAT  | Clust-9  |


### PTK7 aptamers prediction
Through the previous classification of the aptamer family, a total of 2,395 aptamer sequences binding to the PTK7 protein were identified. Based on these sequences, we trained the FCNARRB model（https://github.com/turningpoint1988/fcnarbb）, enabling accurate prediction of whether unknown sequences can bind to the PTK7 protein.

```

python ./aptamer_family_analysis/fcna_trainer.py -train_data ./data/input/ptk7_2cls_new.csv -external_data ./data/input/external_data.csv -output grad_folder_1215-2 

```

### PTK7 aptamers generation
Based on the results of Aptamer Family Classification, aptamer sequences binding and non-binding to PTK7 protein were used to train the RaptGen model ([https://github.com/hmdlab/raptgen](https://github.com/hmdlab/raptgen)) to generate new sequences with potential PTK7 protein binding affinity.

```

python ./aptamer_family_analysisrun_aptamer_training.py -ptk7_sample_path ./data/input/clust1_ptk7.csv -negatibe_sample_path ./data/input/other_sequences.csv 

```









## Example
Based on the previous aptamer sequence family grouping information, the aptamer family abundance matrix was generated from the UMI count matrix of the aptamer sequences.For example,we generated 'motit_need_1w' matrix.

```
#Read the results of **Raw Data Preparation**.
mrna_sgrna<-Read10X("./CRISPR_result/filtered_feature_bc_matrix/")
aptamer<-Read10X("./Aptamer_result/")

#mRNA abundance matrix
SUM159 <- CreateSeuratObject(counts = raw_mrna_sgrna$`Gene Expression`[rowSums(raw_mrna_sgrna$`Gene Expression`)>0,])
#sgRNA abundance matrix
SUM159[["sgRNA"]] <- CreateAssayObject(raw_mrna_sgrna$`CRISPR Guide Capture`)
#top_aptamer abundance matrix
SUM159[["aptamer_1w"]] <- CreateAssayObject(aptamer[top_aptamer_sequence,])
#aptamer family abundance matrix
SUM159[["motif1w"]]<-CreateAssayObject(motif_need_1w)

```
### Step1: Quality control and cell classification
  Before performing cell quality control, ensure that you have a Seurat object that includes three essential components: mRNA, aptamer and motif.

```
SUM159<-cell_quality (SUM159,
                      count_threshold = 100,
                      feature_threshold = 200,
                      percent_mt_threshold = 10,
                      assay = "sgRNA",
                      save_path = NULL)
```

<div align="center">
  <img src="picture/sgrna_qc1.png" alt="annotation" width="45%" height="45%" style="display: inline-block;"/>
  <img src="picture/sgrna_qc2.png" alt="annotation" width="45%" height="45%" style="display: inline-block;"/>
</div>


  In this step, you will assign a gRNA identity to each cell and calculate enrichment ratios using cell gRNA counts, which were assessed in Step 1. This process involves setting thresholds to categorize cell gRNA effectively.

```
SUM159<-cell_gRNA_identity(SUM159,
                           assay='sgRNA',
                           min_count = 200,
                           min_ratio = 0.7)
```

### Step2: Aptamer-target interaction prediction
  In this step, you will predict which proteins are likely bound by the aptamer families. This involves calculating the differential matrix based on the median difference between target cells and NC (Control) cells, filtering out low-difference families and confusing targets（If more than half of the families are ranked in the top three, then we believe that they are gRNAs that easily cause confusing differences.）, and using a Gaussian Mixture Model (GMM) to refine the predictions.

```
predict_result<-predict_apt_pro(SUM159,
                                assay = "motif1w",
                                top_n = 20,
                                save_path = NULL)
```

### Step3: Data visualization
  Visualization is crucial for interpreting and presenting the results of your aptamer-protein binding predictions. In this step, you'll generate plots to help understand the data and findings from the previous analysis.

```
visualize_aptamer_difference(predict_result,'Clust-1')

```
<div align="center">
  <img src="picture/PTK7_Clust-1.png" alt="annotation" style="width: 50%; height: auto;"/>
</div>



