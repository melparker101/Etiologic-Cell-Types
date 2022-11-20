#############################################################
## Run CELLEX on In-House Data
## melodyjparker14@gmail.com - August 22
#############################################################

##############################
# 1 - Import packages
##############################
import csv
import gzip
import os
import scipy.io
import pandas as pd
import numpy as np 
import cellex

##############################
# 1 - Source Files
##############################
in_counts <- "counts_f.mtx.gz"
in_counts_alt <- "counts_alt_f.mtx.gz"
in_metadata = "./metadata_adipo_sn_f.csv"
in_metadata_alt = "./metadata_alt_adipo_sn_f.csv"

out_esmu <- "adipo_sn_lindgren_ens.esmu.csv.gz"
out_esmu_alt <- "adipo_sn_alt_lindgren_ens.esmu.csv.gz"

##############################
# 2 - Start code
##############################
# Read in matrix files
mat = scipy.io.mmread(in_counts)
mat_alt = scipy.io.mmread(in_counts_alt)

# Make panda dataframes
matrix = pd.DataFrame.sparse.from_spmatrix(mat)
matrix_alt = pd.DataFrame.sparse.from_spmatrix(mat_alt)

# Read in col names (without /n added on)
# https://stackoverflow.com/questions/3142054/python-add-items-from-txt-file-into-a-list
# with automatically closes the file afterwards
# Leave empty line in
with open('counts_col_names.txt', 'r') as f:
    header = [line.strip() for line in f]

matrix.columns=header
matrix.iloc[0:5,0:2]
# f.close() not needed
# check if file is closed
f.closed

# Add alt columns
# Leave empty line in
with open('counts_alt_col_names.txt', 'r') as f:
    header_alt = [line.strip() for line in f]

matrix_alt.columns=header_alt
matrix_alt.iloc[0:5,0:5]
# f.close() not needed

# Add rows
# Leave empty line in
with open('counts_row_names.txt', 'r') as f:
    row_names = [line.strip() for line in f]

matrix.index=row_names
matrix.iloc[0:5,0:2]

# Add alt rows
# Leave empty line in
with open('counts_alt_row_names.txt', 'r') as f:
    row_names_alt = [line.strip() for line in f]

matrix_alt.index=row_names_alt
matrix_alt.iloc[0:5,0:2]

# Read in metadata csv files
metadata = pd.read_csv(in_metadata, index_col=0)
metadata_alt = pd.read_csv(in_metadata_alt, index_col=0)

# Run CELLEX
eso = cellex.ESObject(data=matrix, annotation=metadata, verbose=True)
eso.compute(verbose=True)
# Convert to ens
cellex.utils.mapping.human_symbol_to_human_ens(eso.results["esmu"], drop_unmapped=True, verbose=True)
# Write esmu file
eso.results["esmu"].to_csv(out_emsu)
del eso

# Run CELLEX on alt data
eso = cellex.ESObject(data=matrix_alt, annotation=metadata_alt, verbose=True)
eso.compute(verbose=True)
# Convert to ens
cellex.utils.mapping.human_symbol_to_human_ens(eso.results["esmu"], drop_unmapped=True, verbose=True)
# Write esmu file
eso.results["esmu"].to_csv(out_esmu_alt)
del eso
