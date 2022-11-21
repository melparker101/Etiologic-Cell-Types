#############################################################
## Merging Single Cell Count Data
## melodyjparker14@gmail.com - Sept 22
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
 
############################## 
# 2 - Source file
##############################
# make a vector contailing all the file names
# data_dir is the file containing all the counts matrices
data_dir <- 'data_dir'
# list.files(data_dir)
counts_file_names <- list.files("data_dir")
n <- length(counts_file_names)

##############################
# 4 - Start code
##############################
# read in all of the counts matrices from each sample
# make a vector "my_count_matrices" containing count matrix names, e.g. count_matrix_1
my_count_matrices = c()

for(i in 1:n) { 
	name <- paste("count_matrix", i, sep="_")
	my_count_matrices[i] <- name
	file <- counts_file_names[i]
	assign(name, readSparseCounts(counts_file_names[i]))
}
# ls()
# dim(count_matrix_1)

# make the count matrices into seurat objects 
# seurat_ob_1 <- CreateSeuratObject(counts = cm_1, project = "adipo",min.cells = 2, min.features = 200)
# make a vector "my_seurat_objects" containing the seurat object names
my_seurat_objects=c()
for (i in 1:n){
	name <- paste("count_matrix", i, sep="_")
	seurat_ob_name <- paste("seurat_ob", i, sep="_")
	my_seurat_objects[i] <- seurat_ob_name_name
	string <- paste(seurat_ob_name_name,"<- CreateSeuratObject(counts = ",name,", project =\"adipo\", min.cells = 2, min.features = 200)")
	eval(parse(text = string))
}
m <- length(my_seurat_objects)

# make a list of annotation names from the file names
# the sample number is always the same length in this case
# the file extention is also always the same length
my_biosample_ids <- c()
my_biosample_ids <- substring(counts_file_names,12,nchar(counts_file_names)-11)

# this is the command we want to achieve:
# adipo.big <- merge(seurat_ob_name_1, y = y_input, add.cell.ids = my_biosample_ids, project = "adipo")

# can't insert my_seurat_objects vector into merge 
# so create a string "merge_adipo" of seurat objects c(...,...)
# NOT including seurat_ob_name_1
y_input <- paste("c(", paste(my_seurat_objects[2:m],collapse = ","),")",sep="")

# make a string with all of this
merge_adipo <- paste("adipo.big <- merge(seurat_ob_name_1, y = c(", paste(my_seurat_objects[2:m],collapse = ","), "), add.cell.ids = my_biosample_ids, project = \"adipo\")",sep="")
# run string as command
adipo.combined <- eval(parse(text = merge_adipo))

# make merged count matrix from seurat object
counts_sc_merged <- as.data.frame(adipo.combined@assays$RNA@counts)
# have a look
counts_sc_merged[1:10, 1:5])

# convert to matrix
counts_sc_matrix <- as.matrix(counts_sc_merged)
mat <- as(counts_sc_matrix, "dgTMatrix")
