#############################################################
## Filtering NAFLD Summary Statistics Data
## melodyjparker14@gmail.com - Sep 22
#############################################################

###################################
# 0 - WORKFLOW AND SETUP
###################################
# The filtering steps are to remove entries that:
# 1. fail HWE (variants in file ukbb_imputed.hardy.fail);
# 2. have high missingness (variants in file ukbb_imputed.vmiss.fail5;
# 3. markers with minor allele frequency <1% (get this from column AF_Allele2: MAF = AF_Allele2 if this is < 0.5, or 1 - AF_Allele2 if AF_Allele2 ≥ 0.5) (Or you can use AF_Allele2 directly by filtering out variants that have AF_Allele2< 0.01 or AF_Allele2 > 0.99)
# 4. markers with SE> 10
# 5. markers with imputation info score < 0.8 (use column imputationInfo)
# 6. exclude duplicate variants 
# 7. exclude multi-allelic variants

# Results were filtered to get rid of any entries without a rsID
# grep rs UKBB.NAFLD.GWAS.results.txt > UKBB_NALFD.txt

###################################
# 1 - LOAD LIBRARIES
###################################
library(data.table)
library(dplyr)

###################################
# 2 - SOURCE FILES
###################################
input_file <- "UKBB.IMPUTED500k.EUR.NAFLD.ALLchr.SAIGE.GWAS.results.txt"
output_file <- "NAFLD_premunge.txt"
in_HWE <- "ukbb_imputed.hardy.fail"
in_miss <- "ukbb_imputed.vmiss.fail5"

###################################
# 3 - START MAIN CODE
###################################
# Read in data
# df <- read.table(input_file, header = TRUE, sep = "", dec = ".")
# fail_HWE <- read.table(in_HWE, header = TRUE, sep = "", dec = ".")
# fail_miss <- read.table(in_miss, header = TRUE, sep = "", dec = ".")

dt <- fread("UKBB_NALFD.txt")
fail_HWE <- fread(in_HWE,header=FALSE)
fail_miss <- fread(in_miss)

### 1 ###
# make a vector from HWE data
fail_HWE_vec <- as.vector(fail_HWE$V1)
dt2 <- subset(dt, !(rsid %in% fail_HWE_vec))

# dim(dt)
# dim(dt2)

### 2 ###
fail_miss_vec <- as.vector(fail_miss$ID)
dt3 <- subset(dt2, !(rsid %in% fail_miss_vec))

# dim(dt3)

### 3 ###
library(dplyr)
dt4 <- filter(dt3, AF_Allele2 >= 0.01 & AF_Allele2 <= 0.99)

# dim(dt4)

### 4 ###
dt5 <- filter(dt4, SE <= 10)
# dim(dt5)

### 5 ###
dt6 <- filter(dt5, imputationInfo >=0.8)
# dim(dt6)

### 6 ###
dt7 <- dt6[nchar(dt6$Allele1) == 1,]
dt8 <- dt7[nchar(dt7$Allele2) == 1,]

# dim(dt7)
# dim(dt8)

### 7 ###
# Look for dupicates
dt8 %>% group_by(rsid) %>% filter(n()>1) %>% summarize(n=n())
# Munging should get rid of the duplicates

# Write to file
write.table(dt8,output_file,col.name = TRUE, row.names=FALSE,sep="\t", quote = FALSE)
