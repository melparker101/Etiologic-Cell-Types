#############################################################
## Converting counts to csv (only works for small matrices)
## melodyjparker14@gmail.com - July 22
## Based on code from Saskia Reibe
## This code converts a count matrix to a csv file in R
#############################################################

# install.packages("Seurat")
library(dplyr)
library(Seurat)
library(patchwork)

# setwd("J:/Adiposites/CELLEX")
data_dir <- 'matrix_files'

# list.files(data_dir)  # should show 'barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz'
data <- Read10X(data.dir = data_dir)

# filter
alldata = CreateSeuratObject(counts = data, project = "adipo",min.cells = 3, min.features = 200)

# alldata [1:10, 1:3]
# dim(alldata )
# as.data.frame(alldata@assays$RNA@counts[1:10, 1:5])

counts <- GetAssayData(object = alldata, slot = "counts")

# write csv to file
write.csv(counts,"counts.csv")
