

##==========================Input========================##
# path_fpkm='/home/huangwb8/Test/RNA-Seq_2/output/fpkm_cufflinks'

## Note:
# The content of bed file is like:
# chr1 3216021 3216967 Xkr4
# chr1 3421701 3421900 Xkr4
# chr1 3670551 3671347 Xkr4

##==============assistant base function=================##
source("/home/huangwb8/bin/ShellBase.R",encoding = "UTF-8")

getOneGeneFPKM <- function(x){
  if(F){
    x <- df.i2[grep("ENSG00000070614",df.i2$gene_id),]
  }
  
  
  
  
}




##====================environment=======================##
nd.pac=c("dplyr","plyr","readr")
Plus.library(nd.pac)

# grobal options:
options(scipen = 1)

##===========================Programe=======================##

file <- list.files(path = path_fpkm,pattern = "genes.fpkm_tracking",full.names = T,recursive = T) # file=file[1:3]

df <- NULL
for(i in 1:length(file)){ # i=1
  
  f.i <- file[i]
  s <- Fastextra(f.i,"/")
  sample <- s[match("genes.fpkm_tracking",s)-1]
  
  ## get FPKM data
  df.i <- read_delim(f.i,delim = "\t") # colnames(df.i)
  df.i2 <- df.i[c("gene_id","length","FPKM")]
  
  
  colnames(df.i2)[2] <- sample
  
  ## merge data
  if(is.null(df)){
    df <- df.i2
  } else {
    df <- left_join(df,df.i2,by="gene_id")
  }
}
df <- as.data.frame(df,stringsAsFactors = F)

table(duplicated(df$gene_id))


rownames(df) <- Fastextra(df$gene_id,"[.]",1)

df <- as.matrix(df)
rownames(df) <- Fastextra(df[,"gene_id"],"[.]",1)
mt <- matrix(as.numeric(df),nrow = nrow(df),byrow = T,dimnames = list(rownames(df),colnames(df)))


str(df)














