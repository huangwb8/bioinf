

##==========================Input========================##
# path_input = "~/Test/annovar/cnv/input/filtered_DCIS_P10_Exon.vcf.gz"

## Note:
# A gatk filter result: *.vcf.gz

##==============assistant base function=================##
sc <- "/home/huangwb8/bin/ShellBase.R"
source(sc, encoding = "UTF-8")
LuckyVerbose("Load source functions from: ",sc)

##====================environment=======================##
nd.pac=c("dplyr","plyr","readr")
scipen = 1

LuckyVerbose("Load grobal options: scipen=",scipen,", needed package: ",paste0(nd.pac,collapse = ","))
Plus.library(nd.pac)
options(scipen = scipen)

##===========================Programe=======================##

vcf <- read.table(path_input,sep = "\t")
n <- basename(path_input)
n2 <- Fastextra(n,".vcf",1)
colnames(vcf) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",n2)
vcf$ID <- vcf$POS
















