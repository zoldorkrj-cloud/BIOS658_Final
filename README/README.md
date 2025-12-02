# RNA-seq Analysis Workflow for BIOS658 Final (GSE72536)


---
## How to get started

**Clone the repository**

From your terminal or GitHub Desktop: git clone https://github.com/zoldorkrj-cloud/BIOS658_Final.git

**Open the RStudio project**

Navigate into the cloned folder BIOS658_Final/.
Open the file BIOS658_Final.Rproj in RStudio. This sets the working directory to the repo root so all relative paths resolve correctly.

## 1. `01_quality_assessment.R`

**Purpose:**  
Perform quality control checks on raw count data.

**Inputs:**  
- `mtx_clean.tsv` (gene counts matrix with genes as rows, samples as columns)  
- Sample metadata (`gds_df.csv` with GEO accessions and HPV status)

**Outputs:**  
- `qc_library_sizes.png` – barplot of library sizes  
- `qc_density.png` – density plots of log2 counts per sample  
- `qc_boxplot.png` – boxplots of log2 counts distributions  
- `qc_correlation_heatmap.png` – sample correlation heatmap  
- `qc_pca_before.png` – PCA plot before batch correction  

---

## 2. `02_preprocessing.R`

**Purpose:**  
Filter lowly expressed genes and normalize counts using edgeR.

**Inputs:**  
- `mtx_clean.tsv` (raw counts matrix)

**Outputs:**  
- Filtered `DGEList` object (`dge_filtered.rds`)  
- Normalized counts (`log2cpm.tsv`)  

---

## 3. `03_batch_correction.R`

**Purpose:**  
Correct batch effects using SVA.
 
**Inputs:**  
- `log2cpm.tsv` (normalized expression matrix)  
- `meta.csv` (sample metadata with HPV status and optional batch info)

**Outputs:**  
- `pca_before_batch.png` – PCA before correction  
- `pca_after_batch.png` – PCA after correction  
- `hc_before_correction.png` – hierarchical clustering before correction  
- `hc_after_correction.png` – hierarchical clustering after correction  
- `log2cpm_corrected.tsv` – batch-corrected expression matrix  

---

## 4. `04_differential_expression.R`

**Purpose:**  
Perform differential expression analysis using edgeR.

**Inputs:**  
- Filtered and normalized `DGEList` object (`dge_filtered.rds`)  
- `meta.csv` (sample metadata with HPV status)

**Outputs:**  
- `edgeR_DE_HPVs_vs_HPVn.csv` – full DE results table  
- `edgeR_top10_boxplots.png` – boxplots of top 10 DE genes  
- `edgeR_volcano_symbol.png` – volcano plot with gene symbols labeled  

---

## 5. `05_functional_enrichment.R`

**Purpose:**  
Perform functional enrichment analysis (ORA and GSEA) using GO and KEGG.

**Inputs:**  
- `edgeR_DE_HPVs_vs_HPVn.csv` (DE results table with logFC and FDR)  
- Gene annotation (`org.Hs.eg.db`)

**Outputs:**  
- `ORA_GO_BP.png` – GO Biological Process enrichment  
- `ORA_GO_MF.png` – GO Molecular Function enrichment  
- `ORA_GO_CC.png` – GO Cellular Component enrichment  
- `ORA_KEGG.png` – KEGG pathway enrichment  
- `GSEA_fgsea_KEGG.csv` – GSEA results for KEGG pathways  

---

# Notes
- All plots are saved in `data/` by default.   

## Original data set

**GSE72536 – RNA-seq of HPV-positive and HPV-negative head and neck cancers.**
