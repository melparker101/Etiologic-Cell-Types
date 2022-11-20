# CELLEX Workflow for Single-Nuclei Data

### 1. Create a new conda environment and install CELLEX
https://github.com/perslab/CELLEX

### 2. Download metadata and matrix files (barcodes, features, mtx) for both mouse and human cells
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

### 3. Filter metadata 
Use filter_metadata.R to filter metadata and create two separate files for mouse and human.

### 4. Run CELLEX
Use run_CELLEX.py to make count matrix in python and run CELLEX.
The count matrix can also be made in R with Seurat and uploaded to python, however, saving it as a CSV file may be difficult due to how large and sparse the matrix is.
