

##=========================Description====================##
# fast way to get small RNA counts from SPAR result

##==========================Input========================##
#peakSampleID="/home/huangwb8/Project/RIP-Seq_hnRNPA2B1/input/case/peakSampleID.txt"
#path_bam="/home/huangwb8/Project/RIP-Seq_hnRNPA2B1/output/align/hisat2"
#control_symbol="null"

args <- commandArgs(T)
peakSampleID <- args[1]
path_bam <- args[2]
control_symbol <- args[3]


## content of peakSampleID is like:
# name control/treat..
# repeat1	SRR4175308_null	SRR4175306_hnrnpa2b1	SRR4175307_hnrnpa2b1
# repeat2	SRR4175311_null	SRR4175309_hnrnpa2b1	SRR4175310_hnrnpa2b1

##=========================Usage========================##
# shell:
# Rscript macs2_expression.R ${peakSampleID} ${path_bam} ${control_symbol}



##==============assistant base function=================##
sc <- "/home/huangwb8/bin/ShellBase.R"
source(sc, encoding = "UTF-8")

##====================environment=======================##
nd.pac="tidyr"; Plus.library(nd.pac)

##===========================Programe=======================##
psi <- read.table(peakSampleID,header = F,sep = "\t",check.names = F)
e <- NULL
for(i in 1:nrow(psi)){ # i=1
  psi.i <- as.character(as.matrix( psi[i,]))
  name.i <- psi.i[1]; psi.i <- psi.i[2:length(psi.i)]
  path_psi.i <- list.files(path_bam,pattern = ".bam$",full.names = T) %>% .[Fastgrep(psi.i,.)]
  preffix.i <- ifelse(grepl(control_symbol,psi.i),"-c ","-t ")
  e.i <- paste(preffix.i, path_psi.i,sep = "") %>% paste0(.,collapse = " ")
  df.i <- data.frame(names=name.i,expr=e.i)
  e <- rbind(e,df.i)
}
path_e <- gsub(basename(peakSampleID),"",peakSampleID)
write.table(e,paste0(path_e,"macs2_expression"),quote = F,sep = "\t",col.names = F,row.names = F)
LuckyVerbose("macs2_expression done!")

