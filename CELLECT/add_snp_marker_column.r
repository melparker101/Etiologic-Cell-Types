# Code from Laura Wittemans
# df is summary statistic data
# Columns in df include: pos, chr, Allele1, Allele2
n <- nrow(df)
for (i in 1:n){
	a1 <- toupper(df$Allele1[i])
	a2 <- toupper(df$Allele2[i])
	chrom <- df$chr[i]
	position <- df$pos[i]
	alleles <- paste(sort(c(a1,a2)), collapse="_")
	df[i,"key"]<- paste(chrom,":",position,"_",alleles,sep="")
}
