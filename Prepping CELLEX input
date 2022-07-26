# import packages
import csv
import gzip
import os
import scipy.io
import pandas as pd
import numpy as np 
import cellex

# define MEX directory
matrix_dir = "matrix_files"
# read in MEX format matrix as table
mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx.gz"))
 
# list of transcript ids, e.g. 'ENSG00000243485'
features_path = os.path.join(matrix_dir, "features.tsv.gz")
feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
 
# list of gene names, e.g. 'MIR1302-2HG'
# gene_names = [row[1] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
 
# list of feature_types, e.g. 'Gene Expression'
# feature_types = [row[2] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode="rt"), delimiter="\t")]
 
# transform table to pandas dataframe and label rows and columns
matrix = pd.DataFrame.sparse.from_spmatrix(mat)
matrix.columns = barcodes
matrix.insert(loc=0, column="feature_id", value=feature_ids)

matrix.iloc[0:10,0:10]

# cut off the decimals in feature_id
matrix['feature_id'] = matrix['feature_id'].apply(lambda x: x.rsplit('.',1)[0])

# read in csv
metadata = pd.read_csv("./cell_metadata.csv", index_col=0)
metadata.iloc[0:10,0:10]

# matrix.set_index('feature_id')
matrix = matrix.set_index('feature_id')

# run cellex
eso = cellex.ESObject(data=matrix, annotation=metadata, verbose=True)
eso.compute(verbose=True)
