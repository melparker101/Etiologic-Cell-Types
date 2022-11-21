#############################################################
## Converting Chromosome Cordinates to RSIDs
## melodyjparker14@gmail.com - August 22
## This script uses VEP to convert chromosome coordinates to RSIDs.
## The data were actually too large, so this script was not actually used for CELLECT preperation.
## https://github.com/Ensembl/ensembl-vep
#############################################################

# Download VEP files
git clone https://github.com/Ensembl/ensembl-vep
git checkout release/103.1

module load VEP/103.1-GCC-10.2.0
# This will automatically load Perl/5.32.0-GCCcore-10.2.0
# If installing manually, load Perl first (find an acceptable version using 'module avail')

# Make sure the .vep folder is in the right folder - not the home dir
# Move .vep file from home dir to ensemble-vep
cd ~
mv .vep software/ensembl-vep/.vep
# Create a syslink to redirect for storing the cache in .vep 
# as it will search in home dir
ln -s software/ensembl-vep/.vep .vep
# Make a shortcut
VEP=software/ensembl-vep/vep

# Install cache
cd software/ensembl-vep
perl INSTALL.pl
# Say no to updating version
# Say no to installing API
# Choose 457 : homo_sapiens_vep_103_GRCh38.tar.gz

#############################################################
# Copy sumstats directory over
cp -R ~/sumstats/ software/ensembl-vep/

# Make a new VEP input file with the chromosome positions in the right format
# e.g. '1 10511 . A G'
awk '{ print $2, $3, ".", $4, $5 }' PCOS_UKBBall.txt > PCOS_input.txt
# Remove the header
sed -i '1d' PCOS_input.txt

# Run VEP
vep ./vep -i PCOS_input.txt -o PCOS_results.txt \
  --dir_cache .vep \
  --cache  \
  --species "homo_sapiens" \
  --check_existing \
  --vcf \
 --symbol --fields Uploaded_variation,Existing_variation\
  --force_overwrite \
  --no_stats \
  --offline
