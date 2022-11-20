#############################################################
## Running CELLEX
## melodyjparker14@gmail.com - July 22
## This script makes a counts matrix data frame from the three matrix files, 
## reads in the metadata as a data frame and runs CELLEX.
## Default is mouse. Adjust for human cells.
## Tutorial used:
## https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices
#############################################################

# Import packages
import csv
import gzip
import os
import scipy.io
import pandas as pd
import numpy as np 
import cellex

# Define MEX directory and metadata file
matrix_dir = "matrix_files"
metadata_file = "./metadata_mm.csv"

# Read in MEX format matrix as table
mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx.gz"))
 
# List of transcript ids
features_path = os.path.join(matrix_dir, "features.tsv.gz")
feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]

barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode="rt"), delimiter="\t")]
 
# Transform table to pandas dataframe and label rows and columns
matrix = pd.DataFrame.sparse.from_spmatrix(mat)
matrix.columns = barcodes
matrix.insert(loc=0, column="feature_id", value=feature_ids)

# Edit the feature ids
matrix['feature_id'] = matrix['feature_id'].apply(lambda x: x.rsplit('.',1)[0])
# Set index of matrix
matrix = matrix.set_index('feature_id')
matrix.iloc[0:10,0:10]

# Read in metadata
metadata = pd.read_csv(metadata_file, index_col=0)
metadata = metadata.set_index('cell_id')
metadata.iloc[0:10,0:10]

# Run CELLEX
eso = cellex.ESObject(data=matrix, annotation=metadata, verbose=True)
eso.compute(verbose=True)

cellex.utils.mapping.mouse_ens_to_human_ens(eso.results["esmu"], drop_unmapped=True, verbose=True)

# eso.results["esmu"].to_csv("mydataset.esmu.csv.gz")
eso.save_as_csv(keys=["all"], verbose=True)
