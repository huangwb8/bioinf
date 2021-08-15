

##==========================Input========================##
# path_fpkm='/home/huangwb8/Project/RNA-Seq_A2B1-KO/output/fpkm_cufflinks'

## Note:
# The content of bed file is like:
# chr1 3216021 3216967 Xkr4
# chr1 3421701 3421900 Xkr4
# chr1 3670551 3671347 Xkr4

##==============assistant base function=================##
source("/home/huangwb8/bin/ShellBase.R",encoding = "UTF-8")

##====================environment=======================##
nd.pac=c("dplyr","plyr","readr")
Plus.library(nd.pac)

# grobal options:
options(scipen = 1)

##===========================Programe====================##

file <- list.files(path = path_fpkm,pattern = "genes.fpkm_tracking",full.names = T,recursive = T) # file=file[1:3]

df <- NULL
for(i in 1:length(file)){ # i=1
  
  f.i <- file[i]
  s <- Fastextra(f.i,"/")
  sample <- s[match("genes.fpkm_tracking",s)-1]
  
  ## get FPKM data
  df.i <- read_delim(f.i,delim = "\t")
  s.i <- Fastextra(df.i$gene_id, '[.]', 1)
  df.i2 <- mergeMatrixDup(
    as.matrix(df.i[c("FPKM", "FPKM_conf_lo", "FPKM_conf_hi")]), 
    mergeCol = F,
    mergeRow = T,
    refRow = s.i, 
    fun_row = sum)
  df.i3 <- cbind(gene_id = rownames(df.i2), as.data.frame(df.i2))
  colnames(df.i3)[2:4] <- paste(sample, colnames(df.i3)[2:4], sep = '_')
  
  ## merge data
  if(is.null(df)){
    df <- df.i3
  } else {
    df <- left_join(df, df.i3, by="gene_id")
  }
  
}
# df$gene_id[duplicated(df$gene_id)]

df2 <- as.matrix(df[,-1])
rownames(df2) <- as.character(df[,1])
df3 <- df2[,grepl('_FPKM$',colnames(df2))]

saveRDS(df3,file=paste0(path_fpkm,"/matrix_fpkm.rds"))
write.table(df3,paste0(path_fpkm,"/matrix_fpkm.txt"),sep = "\t",col.names = T,row.names = T,quote = F)

LuckyVerbose("All done!")

