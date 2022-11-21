#############################################################
## Merge Counts - Single Cell Atlas Data
## melodyjparker14@gmail.com - August 22
## This script merges count data from 36 individual samples
## Tutorial for merging two samples:
## https://satijalab.org/seurat/articles/merge_vignette.html
#############################################################

##############################
# 1 - Load librairies
##############################
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratData)
library(scuttle)  # for reading in sparse matrices
library(Matrix)
library(data.table)

############################## 
# 2 - Source files
##############################
data_dir <- 'data_dir'
in_metadata <- "GSE176171_cell_metadata.tsv"

out_counts <- "matrix_sc_f.mtx.gz"
out_counts_col_names <- "matrix_f_col_names.txt"
out_counts_row_names <- "matrix_f_row_names.txt"

##############################
# 3 - Define functions
##############################
# Funtion: https://slowkow.com/notes/sparse-matrix/
# Enables saving matrix as a compressed file
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
# 4 - Start main code
##############################
# Make a vector contailing all the file names
# data_dir is the file containing all the counts matrices
list.files(data_dir)
c_matrices <- list.files("data_dir")
n <- length(c_matrices)

# Read in all of the counts matrices from each sample
# Make a vector x containing count matrix names, e.g. cm_1
x=c()
for(i in 1:n) { 
	name <- paste("cm", i, sep="_")
	x[i] <- name
	file <- c_matrices[i]
	assign(name, readSparseCounts(c_matrices[i]))
}
# ls()
# dim(cm_1)

# Make the count matrices into seurat objects 
# so_1 <- CreateSeuratObject(counts = cm_1, project = "adipo",min.cells = 2, min.features = 200)
# Make a vector y containing the seurat object names
y=c()
for (i in 1:n){
	name <- paste("cm", i, sep="_")
	so_name <- paste("so", i, sep="_")
	y[i] <- so_name
	string <- paste(so_name,"<- CreateSeuratObject(counts = ",name,", project =\"adipo\", min.cells = 2, min.features = 200)")
	eval(parse(text = string))
}

# Make a list of annotation names from the file names
# The sample number is always the same length in this case
z <- c()
z <- substring(c_matrices,12,nchar(c_matrices)-11)

# This is the command we want to achieve:
# adipo.big <- merge(so_1, y = y_input, add.cell.ids = z, project = "adipo")

# Can't insert y vector into merge 
# so create a string "merge_adipo" of seurat objects c(...,...)
# NOT including so_1
y_input <- paste("c(", paste(y[2:length(y)],collapse = ","),")",sep="")

# Make a string with all of this
merge_adipo <- paste("adipo.big <- merge(so_1, y = c(", paste(y[2:length(y)],collapse = ","), "), add.cell.ids = z, project = \"adipo\")",sep="")
# Run string as command
adipo.combined <- eval(parse(text = merge_adipo))

# Make merged count matrix from seurat object
counts_sc_merged <- as.data.frame(adipo.combined@assays$RNA@counts)
# counts_sc_merged[1:10, 1:5])

# Convert to matrix
counts_sc_matrix <- as.matrix(counts_sc_merged)
m <- as(counts_sc_matrix, "dgTMatrix")

##############################
# Read in metadata
metadata <- fread(in_metadata)

# Filter to cell_id and cell_type only
metadata <- metadata[,c(1,57)]
# head(metadata)

# Make a vector of cell ids to use for filtering count matrix
cell_ids <- metadata$cell_id

# Filter the count data columns
# ids.use is the overlap of cell ids between metadata and count data
# Make a new filtered matrix mf
ids.use <- colnames(m)[colnames(m) %in% cell_ids]
mf <- m[,ids.use]

# Filter the metadata too
metadata_filtered <- subset(metadata, cell_id %in% ids.use)

# Save the sparse matrix as a compressed file
# Save counts and metadata files
# write.csv(metadata_filtered,"metadata_sc_f.csv")
writeMMgz(mf, out_counts)

mf_header <- colnames(mf)
mf_rows <- rownames(mf)
write.table(mf_header,out_counts_col_names,row.names=FALSE,col.names=FALSE,sep="\t", quote = FALSE)
write.table(mf_rows,out_counts_row_names,row.names=FALSE,col.names=FALSE,sep="\t", quote = FALSE)
