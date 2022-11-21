# https://www.biostars.org/p/312369/
# PCOS_UKBall_m.txt is the sumstats file where scientific notation has already been reverted.

# Install bedtools
# conda install bedtools

# Make a file with only required cols
# common_all_20180423_b37.vcf is the dbSNP VCF file
# chr pos pos
awk '{ print $2, $3, $3; }' PCOSUKall_m.txt > positions.txt

# Intersect
bedtools intersect -a positions.txt -b common_all_20180423_b37.vcf -wa -wb | awk '{ print $4,$5,$6 }' > map_og.txt

###################################
# Alternative:

# Download the Single Nucleotide Polymorphism Database VCF file and convert to BED file
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/
# More recent version (version 155, .25 for b37): https://ftp.ncbi.nlm.nih.gov/snp/latest_release/
# https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
wget -qO- ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz | gunzip -c - | vcf2bed --sort-tmpdir=${PWD} --max-mem=2G - > hg19.dbSNP150.bed

# Make a positions input file (.bed)
# chr pos-1 pos
awk -v OFS="\t" '{ print $2, ($3 - 1), $3; }' PCOS_UKBall_m.txt | sort-bed - > positions_PCOS.bed

# make map file
bedmap --echo --echo-map-id --delim '\t' positions_PCOS.bed hg19.dbSNP150.bed > answer_PCOS.bed
