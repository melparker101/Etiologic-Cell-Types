# CELLEX Workflow
Create a new conda environment and install CELLEX
- https://github.com/perslab/CELLEX

## Atlas Single-Nuclei Data

### 1. Download metadata and matrix files (barcodes, features, mtx) for both mouse and human cells
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176171
```
mkdir mouse_matrix_files
cd mouse_matrix_files

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176171/suppl/GSE176171_Mm10X.counts.barcodes.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176171/suppl/GSE176171_Mm10X.counts.features.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176171/suppl/GSE176171_Mm10X.counts.mtx.gz
cd ..

mkdir human_matrix_files
cd human_matrix_files
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176171/suppl/GSE176171_Hs10X.counts.barcodes.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176171/suppl/GSE176171_Hs10X.counts.features.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176171/suppl/GSE176171_Hs10X.counts.mtx.gz
cd ..

# GSE176171_cell_metadata.tsv.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176171/suppl/GSE176171_cell_metadata.tsv.gz
```

### 2. Filter metadata 
Use atlas_filter_metadata.R to filter metadata and create two separate files for mouse and human.

### 3. Run CELLEX
Use atlas_run_CELLEX.py to make count matrix in python and run CELLEX. Do this separately for mouse and human data.
The count matrix can also be made in R with Seurat and uploaded to python, however, saving it as a CSV file may be difficult due to how large and sparse the matrix is.

## In-House Data
The data were annotated in two different ways, thus two metadata files were produced.

### 1. Filter count data for each metadata file
Use the lindg_adipo_data_prep.R script.

### 2. Run CELLEX
Read in the count data and metadata to python and then run CELLEX. Use lindg_run_CELLECT.py script.
