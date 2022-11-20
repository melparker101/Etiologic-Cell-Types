#############################################################
## Preprocessing Metadata
## melodyjparker14@gmail.com - July 22
## This code filters the atlas metadata into CELLECT-input format, separating for mouse and human cells.
## The metadata contains information for both mouse and human cells.
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176171
#############################################################

# Read in metadata and barcodes
metadata <- read.csv("metadata.csv")
barcodes_mm <- read.table("barcodes_mm.tsv")
barcodes_mm <- barcodes_mm[,1]
barcodes_hs <- read.table("barcodes_hs.tsv")
barcodes_hs <- barcodes_hs[,1]

# colnames(metadata)
# head(metadata[,c(1,55:58)])

# Filter the metadata to only contain relevant columns (to fit the input shape for CELLEX)
metadata <- metadata[,c(1,57)]
colnames(metadata) <- c("cell_id","cell_type")

# Filter to only include cells in each barcode list
metadata_mm <- metadata[match(barcodes_mm,metadata$cell_id),]
metadata_mm <- metadata[match(barcodes_mm,metadata$cell_id),]

# Could also use this
# metadata_mm <- metadata[grep("mm", metadata$cell_id),]
# metadata_hs <- metadata[grep("hs", metadata$cell_id),]

# Check the dimensions match
length(barcodes_mm)
dim(metadata_mm)
length(barcodes_hs)
dim(metadata_hs)

# Write to csv
write.csv(metadata_mm, "metadata_mm.csv")
write.csv(metadata_hs, "metadata_hs.csv")
