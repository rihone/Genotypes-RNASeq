library(vcfR)
#library(adegenet)
#library(adegraphics)
#library(pegas)
#library(StAMPP)
#library(lattice)
#library(gplots)
#library(ape)
#library(ggmap)

## PART1 ##
setwd("/Volumes/michoel_grp/HipSci/Genotype/Open")
filenames <- list.files(path = ".", pattern = "vcf.gz$", recursive = TRUE)

#setwd("~/Downloads/EGA_download_client_2/request_EGAD00010001328")
#filenames <- Sys.glob("*.vcf.gz")

# VCF to Genlight
for (i in 1:length(filenames)) {
  currentfile <- filenames[i]
  outfile <- gsub(".vcf.gz", ".txt", currentfile)
  vcf <- read.vcfR(currentfile)
  aa.genlight <- vcfR2genlight(vcf, n.cores=8)
  snpmatrix_out <- as.matrix(aa.genlight)
  snpmatrix_out <- t(snpmatrix_out)
  write.table(snpmatrix_out, file = outfile, col.names=NA)
}

# VCF -> SNP_ID/CHR/POS
for (i in 1:length(filenames)) {
  currentfile <- filenames[i]
  outfile <- gsub(".vcf.gz", ".snpids", currentfile)
  vcf <- read.vcfR(currentfile)
  vcfattribute <- vcf@fix[,1:5]
  snpmatrix_out <- as.matrix(vcfattribute)
  write.table(snpmatrix_out, file = outfile, quote = FALSE, row.names = FALSE)
}


## PART 2 ##
setwd("/Volumes/michoel_grp/HipSci/Processed_data/Genotypes")
filenames <- list.files(path = ".", pattern = ".snpids$", recursive = TRUE)
# Concatenate SNP Attributes
library(dplyr)
currentfile <- filenames[1]
snps <- read.table(currentfile, header = TRUE)
snps <-snps %>% as.data.frame() %>% mutate(UNI = paste(!!!rlang::syms(c("CHROM", "POS", "REF", "ALT")), sep="_"))
for (i in 2:length(filenames)) {
  currentfile <- filenames[i]
  newfile<-read.table(currentfile, header = TRUE)
  newfile <- newfile %>% as.data.frame() %>% mutate(UNI = paste(!!!rlang::syms(c("CHROM", "POS", "REF", "ALT")), sep="_"))
  snps <- rbind(snps, newfile)
  snps <- snps[!duplicated(snps$UNI),]
}
snpmatrix_out <- as.matrix(snps[,c(3,1,2,4,5)])
write.table(snpmatrix_out, file = "all_snps_attribute_open2.txt", quote = FALSE, row.names = FALSE)

#unix
sort -k2,2n -k3,3 match_snps_attribute_open.txt > match_snps_attribute_open2.txt 

curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp150.txt.gz" | gunzip -c | cut -d$'\t' -f 2,3,4,5 | head





