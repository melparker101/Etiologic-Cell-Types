#############################################################
## In-House Adipose Data CELLEX preperation
## melodyjparker14@gmail.com - Sept 22
## This script filters the count data for two different annotation files.
## The count data is then saved as a matrix file with it's column and row names saved separately.
## There must be a better way to do this, but it works.
#############################################################

##############################
# 1 - Load librairies
##############################
library(data.table)
library(Matrix)

############################## 
# 2 - Source files
##############################
in_counts_file <- "adipo_sc_counts_no_normal.txt"
in_metadata_file <- "meta_data_adipo_sn.txt"
in_metadata_alt_file <- "alt_meta_adipo_sn.txt"

out_counts_file <- "counts_f.mtx.gz"
out_counts_alt_file <- "counts_alt_f.mtx.gz"
out_metadata <- "metadata_adipo_sn_f.csv"
out_metadata_alt <- "metadata_alt_adipo_sn_f.csv"
out_counts_col_names <- "counts_col_names.txt"
out_counts_row_names <- "counts_row_names.txt"
out_counts_alt_col_names <- "counts_alt_col_names.txt"
out_counts_alt_row_names <- "counts_alt_row_names.txt"

##############################
# 3 - Define functions
##############################
# Funtion: https://slowkow.com/notes/sparse-matrix/
writeMMgz <- function(x, file) {
  mtype <- "real"
  if (is(x, "ngCMatrix")) {
    mtype <- "integer"
  }
  writeLines(
    c(
      sprintf("%%%%MatrixMarket matrix coordinate %s general", mtype),
      sprintf("%s %s %s", x@Dim[1], x@Dim[2], length(x@x))
    ),
    gzfile(file)
  )
  data.table::fwrite(
    x = summary(x),
    file = file,
    append = TRUE,
    sep = " ",
    row.names = FALSE,
    col.names = FALSE
  )
}

##############################
# 3 - Start main code
##############################

# Read in data
counts <- as.data.frame(fread(counts_file))
metadata <- fread(metadata_file)
metadata_alt <- fread(metadata_alt_file)

# This changes - to .
# counts_df <- as.data.frame(counts)
# counts <- data.frame(counts[,-1], row.names = counts[,1])

# Make first column of counts row names
row.names(counts) <- counts$V1
counts <- counts[-1]

# Filter count data
cell_ids <- metadata$cell_id
ids.use <- colnames(counts)[colnames(counts) %in% cell_ids]
counts_filtered <- counts[,ids.use]

cell_ids_alt <- metadata_alt$cell_id
ids_alt.use <- colnames(counts)[colnames(counts) %in% cell_ids_alt]
counts_alt_filtered <- counts_df[,ids_alt.use]

# Filter metadata
metadata_filtered <- subset(metadata, cell_id %in% ids.use)
metadata_alt_filtered <- subset(metadata_alt, cell_id %in% ids_alt.use)

# Convert count data to dgTMatrix
counts_f_matrix <- as.matrix(counts_filtered)
counts_alt_f_matrix <- as.matrix(counts_alt_filtered)
m <- as(counts_f_matrix, "dgTMatrix")
m_a <- as(counts_alt_f_matrix, "dgTMatrix")

# Write matrix
writeMMgz(m, out_counts_file)
writeMMgz(m_a, out_counts_alt_file)

# Write metadata
write.csv(metadata_filtered, out_metadata)
write.csv(metadata_filtered, out_metadata_alt)

# Extract column and row names
counts_col_names <- colnames(counts_filtered)
counts_row_names <- rownames(counts_filtered)
counts_alt_col_names <- colnames(counts_alt_filtered)
counts_alt_row_names <- rownames(counts_alt_filtered)

# Write col and row names to TXT files
write.table(counts_col_names,out_counts_col_names,row.names=FALSE,col.names=FALSE,sep="\t", quote = FALSE)
write.table(counts_row_names,out_counts_row_names,row.names=FALSE,sep="\t",col.names=FALSE, quote = FALSE)
write.table(counts_alt_col_names,out_counts_alt_col_names,row.names=FALSE,col.names=FALSE,sep="\t", quote = FALSE)
write.table(counts_alt_row_names,out_counts_alt_row_names,row.names=FALSE,col.names=FALSE,sep="\t", quote = FALSE)
