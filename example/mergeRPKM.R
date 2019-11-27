


##==========================Input========================##
# path_rpkm='/home/huangwb8/Test/SCS_2/output/bam2rpkm'
# path_bed='/data/reference/CCDS/mouse/mm10.exon_sub1.bed'

## Note:
# The content of bed file is like:
# chr1 3216021 3216967 Xkr4
# chr1 3421701 3421900 Xkr4
# chr1 3670551 3671347 Xkr4


##==============assistant base function=================##
source("/home/huangwb8/bin/ShellBase.R",encoding = "UTF-8")

# get RPKM for one gene
getOneGeneRPKM <- function(x){
  if(F){
    expr = df_rpkm_2
    gene = "Lypla1"
    x <- expr[expr$gene %in% gene,]
  }
  Li <- as.numeric(x[,3]) - as.numeric(x[,2]) # length = end - start
  L <- sum(Li)
  mt_rpkm <- x[,5:ncol(x)]
  mt_rpkm2 <- (Li * mt_rpkm)/L
  RPKM <- colSums(mt_rpkm2)
  RPKM <- as.data.frame(t(as.matrix(RPKM)),stringsAsFactors=F)
  return(RPKM)
}

# merge with getOneGeneRPKM
getRPKM <- function(df,geneCol="gene"){
  
  if(F){
    df <- df_rpkm_2
  }
  
  ## merge for every gene
  rpkm <- ddply(df,.variables = geneCol,getOneGeneRPKM)
  rownames(rpkm) <- as.character(rpkm[,geneCol])
  rpkm <- rpkm[-match(geneCol,colnames(rpkm))]
  return(rpkm)
}


##====================environment=======================##
nd.pac=c("dplyr","plyr")
Plus.library(nd.pac)

# grobal options:
options(scipen = 1)

##===========================Programe=======================##



# get files
file <- list.files(path = path_rpkm,pattern = "rpkm.txt$",full.names = T) # file <- file[1:10]

df <- NULL
for(i in 1:length(file)){ # i=1
  f.i <- file[i]
  df.i <- read.table(f.i,sep = "\t")
  base <- unlist(strsplit(basename(f.i),".rpkm.txt"))
  colnames(df.i) <- c("position",paste0(base,"_reads"),paste0(base,"_RPKM"))
  if(is.null(df)){
    df <- df.i
  } else {
    df <- left_join(df,df.i,by="position")
  }
}

# get read matrix
s_reads <- colnames(df)[grep("_reads",colnames(df))]
df_reads <- df[c("position",s_reads)]
colnames(df_reads)[2:ncol(df_reads)] <- Fastextra(s_reads,"_reads",1)

s_rpkm <- colnames(df)[grep("_RPKM",colnames(df))]
df_rpkm <- df[c("position",s_rpkm)]
colnames(df_rpkm)[2:ncol(df_rpkm)] <- Fastextra(s_rpkm,"_reads",1)


# bed annotation
bed <- read.table(path_bed,sep = "\t",stringsAsFactors = F)
colnames(bed) <- c("chr","start","end","gene")
bed$position <- as.numeric(rownames(bed))
df_rpkm_2 <- left_join(bed,df_rpkm,by="position")
df_reads_2 <- left_join(bed,df_reads,by="position")
df_rpkm_2 <- df_rpkm_2[-match("position",colnames(df_rpkm_2))]
df_reads_2  <- df_reads_2 [-match("position",colnames(df_reads_2))]

# merge exon to gene——RPKM
RPKM <- getRPKM(df_rpkm_2)
RPKM <- RPKM[factor(rownames(RPKM)),]

# merge exon to gene——Reads Count
geneidfactor <- factor(df_reads_2$gene)
df_reads_3 <- df_reads_2[5:ncol(df_reads_2)]
Counts <- apply(df_reads_3,2,function(x)tapply(x,geneidfactor,sum))

# save matrix
save(RPKM,file = paste0(path_rpkm,"/rpkm_matrix.rda"))
save(Counts,file = paste0(path_rpkm,"/read_matrix.rda"))
write.table(RPKM,paste0(path_rpkm,"/rpkm_matrix.txt"),sep = "\t",col.names = T,row.names = T,quote = F)
write.table(Counts,paste0(path_rpkm,"/read_matrix.txt"),sep = "\t",col.names = T,row.names = T,quote = F)
write.table(df_rpkm_2,paste0(path_rpkm,"/raw_rpkm.txt"),sep = "\t",col.names = T,row.names = F,quote = F)
write.table(df_reads_2,paste0(path_rpkm,"/raw_read.txt"),sep = "\t",col.names = T,row.names = F,quote = F)

