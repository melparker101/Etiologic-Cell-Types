# Example munging script

python ldsc/mtag_munge.py \
--sumstats sumstats/check/bmi_premunge.txt \
--merge-alleles data/ldsc/w_hm3.snplist \
--a1 Tested_Allele \
--a2 Other_Allele \
--keep-pval \
--p P \
--out sumstats/check/BMI_munged
