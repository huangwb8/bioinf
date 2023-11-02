

##=====================Information========================##
# Version: 0.1.0
# Author: Weibin Huang
# Merge samples' stringtie results into an TPM expression matrix and output the matrix as *.txt and *.rds.


##========================Usage===========================##
# Rscript ~/bin/example/mergeStringtieTPM.R --help
# nohup Rscript ~/bin/example/mergeStringtieTPM.R --path_stringtie ${path_res}> ${path_log}/mergeStringtieTPM.log 2>&1 &
 
# tree -L 2 ${path_res} is like: 

# ├── SRR11873614.Aligned.sortedByCoord.out
# │   ├── cov_refs.gtf
# │   ├── e2t.ctab
# │   ├── e_data.ctab
# │   ├── gene_abund.tab
# │   ├── i2t.ctab
# │   ├── i_data.ctab
# │   ├── merged.gtf
# │   └── t_data.ctab
# ├── SRR11873614.Aligned.sortedByCoord.out.transcript.tpm.tab
# ├── SRR11873615.Aligned.sortedByCoord.out
# │   ├── cov_refs.gtf
# │   ├── e2t.ctab
# │   ├── e_data.ctab
# │   ├── gene_abund.tab
# │   ├── i2t.ctab
# │   ├── i_data.ctab
# │   ├── merged.gtf
# │   └── t_data.ctab
# ├── SRR11873615.Aligned.sortedByCoord.out.transcript.tpm.tab
# ...

##==============assistant base function=================##
source("/home/huangwb8/bin/ShellBase.R",encoding = "UTF-8")

##====================environment=======================##
nd.pac=c("dplyr","plyr","readr","optparse")
Plus.library(nd.pac)

# grobal options:
options(scipen = 100)

# Define Shell parameters
# path_stringtie='~/Project/RNA-Seq_A2B1-KO/output/stringtie/hisat2'
option_list <- list(
	make_option(c("-p", "--path_stringtie"), type = "character", default = NULL, help = "The path of stringtie output.")
)

# Shell to R
opt <- parse_args(OptionParser(option_list = option_list))

# R parameters
path_stringtie = opt$path_stringtie

##===========================Programe====================##

f <- list.files(path_stringtie,pattern = 'gene_abund.tab',full.names = T,recursive = T)

df_tpm <- NULL; df_fpkm <- NULL
for(i in 1:length(f)){ # i=1
  
  f.i <- f[i]
  n.i <- rev(Fastextra(f.i,'[/]'))[2]
  df.i <- read.table(f.i,header = T,check.names = F,sep = '\t')
  
  gene.i <- Fastextra(df.i$`Gene ID`, '[.]',1)
  # table(duplicated(gene.i))
  # gene.i[duplicated(gene.i)]
  tpm.i <- tapply(df.i[,'TPM'], gene.i, sum)
  fpkm.i <- tapply(df.i[,'FPKM'], gene.i, sum)
  
  # TPM
  df.i2 <- data.frame(Gene_ID = names(tpm.i), TPM = tpm.i, stringsAsFactors = F,row.names = 1:length(tpm.i))
  colnames(df.i2)[2] <- n.i
  if(is.null(df_tpm)){
    df_tpm <- df.i2
  } else {
    df_tpm <- full_join(df_tpm,df.i2,by = "Gene_ID")
  }
  
  
  # FPKM
  df.i3 <- data.frame(Gene_ID = names(fpkm.i), FPKM = fpkm.i, stringsAsFactors = F,row.names = 1:length(fpkm.i))
  colnames(df.i3)[2] <- n.i
  if(is.null(df_fpkm)){
    df_fpkm <- df.i3
  } else {
    df_fpkm <- full_join(df_fpkm,df.i3,by = "Gene_ID")
  }
  
  
  # TPM for SUPPA
  # df.i2.m <- df.i2[-1]; rownames(df.i2.m) <- df.i2[,1]
  # write.table(df.i2.m, paste0(Fastextra(f.i,n.i,1),n.i,'.tpm.tab'),sep = "\t",col.names = T,row.names = T,quote = F)
  if(T){
    
    x <- readr::read_delim(paste0(path_stringtie,'/',n.i,'/merged.gtf'), delim = '\t', comment = "#",col_names = F) %>% as.data.frame()
    x2 <- x[x[,3] == 'transcript',]
    
    a <- apply(as.matrix(x2[,9]), 1, function(i){
      # i = 'gene_id \"ENSG00000228794.9\"; transcript_id \"ENST00000671208.1\"; ref_gene_name \"LINC01128\"; cov \"0.0\"; FPKM \"0.000000\"; TPM \"0.000000\";'
      i2 <- Fastextra(i,'\"')
      i3 <- i2[seq(2,12,2)]
      # names(i3) <- gsub(';| ','',i2[seq(1,11,2)])
      return(i3)
    })
    a2 <- t(a); 
    colnames(a2) <- c("gene_id","transcript_id" ,"ref_gene_name","cov","FPKM","TPM") # View(a2[1:10,])
    # rownames(a2) <- Fastextra(a2[,"transcript_id"], "[.]",1)
    rownames(a2) <- a2[,"transcript_id"]
    a2 <- as.data.frame(a2)
    a3 <- a2["TPM"]; colnames(a3) <- n.i
    a3[,1] <- as.numeric(as.character(a3[,1]))
    a3[,1][is.na(a3[,1])] <- 0
    # a3 <- na.omit(a3) # table(is.na(a3$D1)); a3[is.na(a3$D1),]

    # Output
    write.table(a3, paste0(Fastextra(f.i,n.i,1),n.i,'.transcript.tpm.tab'),sep = "\t",col.names = T,row.names = T,quote = F)
  }
  
}

# table(duplicated(df_tpm$Gene_ID))
# table(duplicated(df_fpkm$Gene_ID))

# Output
df_tpm2 <- as.matrix(df_tpm[,-1]); rownames(df_tpm2) <- df_tpm[,1]
for(j in 1:ncol(df_tpm2)) df_tpm2[,j][is.na(df_tpm2[,j])] <- 0
saveRDS(df_tpm2, paste0(path_stringtie,'/stringtie.tpm.rds'))
write.table(df_tpm2, paste0(path_stringtie,"/stringtie.tpm.txt"),sep = "\t",col.names = T,row.names = T,quote = F)

df_fpkm2 <- as.matrix(df_fpkm[,-1]); rownames(df_fpkm2) <- df_fpkm[,1]
for(j in 1:ncol(df_fpkm2)) df_fpkm2[,j][is.na(df_fpkm2[,j])] <- 0
saveRDS(df_fpkm2, paste0(path_stringtie,'/stringtie.fpkm.rds'))
write.table(df_fpkm2, paste0(path_stringtie,"/stringtie.fpkm.txt"),sep = "\t",col.names = T,row.names = T,quote = F)

# grobal options:
options(scipen = 1)

