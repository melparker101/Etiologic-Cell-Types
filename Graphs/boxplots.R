#############################################################
## Creating boxplots from CELLECT results
## melodyjparker14@gmail.com - Sept 22
## Plot size needs adjusting for different number of phenotypes.
#############################################################

##############################
# 0 - Load librairies
##############################
library(dplyr)  # for mutate 

############################## 
# 1 - Choose data
##############################
# edit this based on chosen data/phenotypes

# phenotypes need to be consistent with CELLECT output directory names
# for the alternative annotations, add "_ALT" on the end
phenotypes <- c("BMI",  "CAD",  "HT",  "NAFLD",  "PCOS", "T2DadjBMI", "Testosterone_F",  "Testosterone_M",  "WHRadjBMI", "HOMAIRadjBMI")

# data_source options: "atlas", "lindg"
data_source <- "atlas"

# tissue_type options: "adipo", "ovary"
tissue_type <- "adipo"

# RNAseq_type options: "sn", "sc", "mm"
RNAseq_type <- "sn"

############################## 
# 2 - Source file
##############################
# paste general results data file
data_type <- paste("CELLECT", toupper(data_source), toupper(tissue_type), toupper(RNAseq_type), sep="-")
# "CELLECT-ATLAS-ADIPO-SN"

# input files need to be inside loop in main code
output_file <- paste("CELLECT/plots/", data_source, "_", tissue_type, "_", RNAseq_type, sep="")

##############################
# 3 - Define functions
##############################
# function for adding -log(p value) column
results_col <- function(x) {
  mutate(prioritization, "{x}" := -log10(pvalue))
}

##############################
# 4 - Start code
##############################
# make matrix for barplot input
for (pheno in phenotypes){
  # read in results
  prioritization <- read.csv(paste("CELLECT", data_type, pheno, "CELLECT-LDSC/results/prioritization.csv", sep="/"))
  # use function to add -log(pvalue) column
  prioritization <- results_col(pheno)
  # add column to joint df
  if (exists('df1')){
    df1 <- full_join(df1, prioritization[,c(3,7)], by = "annotation")
  } else {
    df1 <- prioritization[,c(3,7)]    
  }
}

df1 <- df1[order(tolower(df1$annotation)),]
rownames(df1) <- df1[,1]
df1 <- df1[,-1]
matrix <- as.matrix(df1)

# choose bonferroni-adjusted significance threshold
n <- nrow(matrix)
m <- ncol(matrix)
sig <-0.05
b_sig = -log10((sig/n))
y_lim <- ceiling(max(matrix))

# generate barplot
pdf(file=output_file, width=20)
{par(lwd = 2)
barplot(matrix, space=c(0,4), width = 1, border=T, beside = TRUE, col = terrain.colors(n), xlim = c(6, y_lim +190), ylim = c(0, y_lim), legend = TRUE, args.legend = list(ncol = 4, bty = "n", xjust=1, x = "top"))
abline(h=b_sig,lty = 2)
title(ylab = '-log(p value)', line=2.2,font.lab=1)
title(main = 'CELLECT - Atlas snRNA-seq Data (Replicate and Other Phenotypes)')
# rect(-100,0,22,3,col = rgb(0,0,1,0), border=F)
# rect(22,0,56,3,col = rgb(0,0,0,1/15), border=F)
abline(h=0)
}
dev.off()
