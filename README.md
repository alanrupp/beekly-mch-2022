# Beekly et al. 2023
Identifying subpopulations of lateral hypothalamic MCH neurons that might use distinct neurotransmitters. To do so, I'm using 3 published scRNA-seq datasets:
* Mickelsen *et al* ([GSE125605](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125065))
* Rossi *et al* ([GSE130597](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130597))
* Jiang *et al* ([GSE146020](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146020))

This analysis isolates MCH neurons from these datasets and integrates them using the Seurat CCA method.

Note: Github file size restrictions prevent uploading all the data files. Please get them from GEO (using the above links) and structure them in a **data** folder like below:

# Files
```
+-- data
    +-- GSE125065
        +-- GSM3562050_AJ17001_barcodes.tsv.gz
        +-- GSM3562050_AJ17001_genes.tsv.gz
        +-- GSM3562050_AJ17001_matrix.mtx.gz
        +-- ...
    +-- GSE130597
        +-- GSM3744476_MLB003_AJP1-N701_merged_GE_clean.bam.dge_1920.txt.gz
        +-- GSM3744477_MLB003_AKP1-N702_merged_GE_clean.bam.dge_2280.txt.gz
        +-- GSM3744478_MLB003_ALP1-N703_merged_GE_clean.bam.dge_2640.txt.gz
        +-- ...
    +-- GSE146020
        +-- GSM4354862_cluster_markers.xlsx
        +-- GSM4354862_count_matrix.csv.gz
        +-- GSM4354862_filtered_feature_bc_matrix.h5
    +-- lha-2019_neuronal.rds
+-- scripts
    +-- figures.R
    +-- GSE125065.R
    +-- GSE130597.R
    +-- GSE146020.R
    +-- integrate.R
    +-- neurotransmitters.R
+-- environment.yml
+-- README.md
```

## Pipeline
This is best run with the `anaconda` package manager to ensure reproducibility. If you don't have `anaconda` installed, you can get it [here](https://www.anaconda.com/products/individual). Or you can use the `environment.yml` file to find the specific packages and versions required to run this analysis.

### `conda` setup
1. Create `conda` environment: `conda env create -f environment.yml`
2. Activate `conda` environment: `conda activate MCH`

### Analysis
1. Isolate MCH neurons from Mickelsen et al. dataset: `Rscript scripts/GSE125065.R`
2. Isolate MCH neurons from Rossi et al. dataset: `Rscript scripts/GSE130597.R`
3. Isolate MCH neurons from Jiant et al. dataset: `Rscript scripts/GSE146020.R`
4. Merge datasets and analyze for neurotransmitters: `Rscript scripts/integrate.R`
5. Find neurotransmitter expression across cells and dataset: `Rscript scripts/neurotransmitters.R`
6. Generate figures for the manuscript: `Rscript scripts/figures.R`

Deactivate `conda` environment: `conda deactivate`

## Results
All analysis results are found in the `/results` directory.
* **GSE125605.rds**: RDS file of the Seurat object containing MCH neurons from GSE125605
* **GSE130597.rds**: RDS file of the Seurat object containing MCH neurons from GSE130597
* **GSE146020.rds**: RDS file of the Seurat object containing MCH neurons from GSE146020
* **integrated.rds**: RDS file of the Seurat object of all datasets integrated

Outputs to generate figures in the manuscript are in the `/figures` directory.
* **figure.svg**: `Inkscape` SVG of the figure
