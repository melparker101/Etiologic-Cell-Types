# Running CELLECT
Workflow for running CELLECT:
- Install condaforge or miniforge and add it to the END of the path in .bashrc
- Clone CELLECT
- Organise file structure
- Prepare ESMU CSV file (CELLEX output)
- Munge GWAS sumstats file

### Install CELLECT
Follow the instructions from the CELLECT github.
https://github.com/perslab/CELLECT < \b>

### Clone CELLECT
```
module use -a /apps/eb/dev/{skylake,ivybridge}/modules/all
module load Anaconda3/2022.05
git lfs env 
git clone --recurse-submodules https://github.com/perslab/CELLECT.git
cd CELLECT
```

### File Structure
CELLECT is the working directory.
```
mkdir EMSU_files  # To contain EMSU input files
mkdir sumstats  # To contain sumstats input files
mkdir config_files  # Containing the config_files
mkdir CELLECT-ATLAS  # Output directory
```

### Generate EMSU files
Use CELLEX. See https://github.com/melparker101/Adipose-Tissue---Etiologic-Cell-Types/tree/main/CELLEX for this workflow.
The CELLEX results must not have any spaces or special characters in the column names, and can only contain one underscore.

### Munge sumstats files
Follow the instructions on the CELLECT github. Make a munging conda environme which installs python2, which is required for the munging script given.
Some tips for running the munge script:
- The SNP column must contain rsIDs. If the file only contains chromosome number and position then these must be converted to rsIDs. If the rsID column is in the format rs140052487:C:A then it must be trimmed (use trim_rsid.sh). This can be done in R using a mapping file generated from a dbSNP VCF (create_map_bedtools.sh) or using VEP (VEP.sh). For large files, VEP is not recommended.
- For dichotomous traits, the effective sample size column Neff must be re-labled as N. The munge script will calculate Neff and label it as N for sumstats that contain N-case and N-controls columns; otherwise, use neff_to_n.sh.
- If sumstats do not contain N, a contant can be added in the munge script using --n-value.

```
# Create munging environment
conda env create -f ldsc/environment_munge_ldsc.yml
conda activate munge_ldsc
```
### Config file
Modify the config.yml file, specifying the paths of:
- input - ESMU CSV file, e.g. EMSU_files/atlas_sc_esmu.csv
- input - GWAS sumstats file, e.g. sumstats/BMI_munged.sumstats.gz
- output - results directory, e.g. CELLECT-ATLAS/BMI 

```
# Example:
BASE_OUTPUT_DIR: CELLECT-ATLAS-BMI

SPECIFICITY_INPUT:
  - id: Hh_esmu
    path: atlas/Hh_esmu.csv

GWAS_SUMSTATS:
  - id: BMI_Pulit2018
    path: atlas/BMI_Pulit2018.sumstats.gz
```

### Run 
Create a conda environment for running snakemake. Run CELLECT via snakemake. The config file will contain all of the input and output information for each run.
```
conda deactivate
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake

snakemake --use-conda -j -s cellect-ldsc.snakefile --configfile config_files/config.yml
```
