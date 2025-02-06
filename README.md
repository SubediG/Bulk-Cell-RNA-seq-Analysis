# Bulk-Cell-RNA-seq-Analysis

This repository contains the code and analysis for exploring differences in gene expression between metastatic **Triple-Negative Breast Cancer (mTNBC) patients and healthy donors**. The project aims to identify differentially expressed genes (DEGs) and uncover molecular insights into mTNBC using transcriptomic data. The analysis is inspired by findings from the article referenced by PMID: 39843922 and utilizes the DESeq2 package in R for differential expression analysis.

**📂 Dataset**

The dataset used in this analysis is publicly available and can be accessed via the GEO database under accession number GSE264108. It includes RNA-seq read counts for 14 samples:
7 samples from metastatic Triple-Negative Breast Cancer (mTNBC) patients.
7 samples from healthy donors.



**🔍 Project Overview**

The goal of this project is to:

i.) Identify differentially expressed genes (DEGs) between mTNBC patients and healthy donors.

ii.) Explore the biological significance of these genes in the context of cancer progression, immune evasion, and metastasis.

iii.) Visualize transcriptional profiles to distinguish mTNBC from healthy samples.

iv.) Lay the groundwork for future machine learning models to predict mTNBC progression.

## 🛠️ Tools & Technologies  

|**Tools**                      | **Description**|
|---------------------------|---------------------------------------------------------------|
| **R Programming**         | Used for data preprocessing, analysis, and visualization. |
| **DESeq2**                | A Bioconductor package for differential gene expression analysis, especially for bulk-cell RNA-seq data. |
| **Volcano Plot** | Helps identify significantly upregulated and downregulated genes. |
| **Heatmap** | Highlights distinct transcriptional profiles between mTNBC and healthy samples. |
| **Git**                   | Used for version control and collaboration. |



**📊 Key Findings**

**1. Differentially Expressed Genes (DEGs)**
Identified a significant number of genes that are either upregulated or downregulated in mTNBC compared to healthy samples.
These genes are potential candidates for further functional studies and could serve as biomarkers for mTNBC.

**2. Top Genes**
Upregulated Genes: Several genes associated with cancer progression, immune evasion, and metastasis were found to be upregulated.
Downregulated Genes: Genes involved in normal cellular functions and tumor suppression were downregulated.

**3. Biological Insights**
A heatmap visualization clearly distinguished between mTNBC and healthy samples, highlighting distinct transcriptional profiles.

This reinforces the idea that mTNBC has a unique molecular signature that could be targeted therapeutically.

**🚀 Next Steps**

The next phase of this project involves leveraging machine learning techniques to predict the presence and progression of mTNBC.
Key steps include:
Incorporating gene expression patterns as features for predictive modeling.

Exploring supervised learning algorithms such as Random Forest and Support Vector Machines (SVM) to classify patients based on gene expression profiles.

Developing models capable of identifying high-risk patients and understanding molecular factors driving tumor progression.


**📄 Usage**
To reproduce the analysis:

```bash
git clone https://github.com/your-username/mTNBC-Transcriptomic-Analysis.git
```

```R
install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("dplyr")
install.packages("pheatmap")
``` 

