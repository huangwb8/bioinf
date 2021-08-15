


##==========================Input========================##
# path_count='/home/huangwb8/Test/RNA-Seq_2/output/count_feature'

## Note:
# The content of bed file is like:
# chr1 3216021 3216967 Xkr4
# chr1 3421701 3421900 Xkr4
# chr1 3670551 3671347 Xkr4

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

## find counts file
file <- list.files(path_count,pattern = "count$",full.names = T)# file <- file[1:3]

## get counts matrix
l <- list();ID <- NULL
for(i in 1:length(file)){ # i=2
  
  f.i <- file[i]
  df.i <- read_delim(f.i,delim = "\t",comment = "#")
  sample <- Fastextra(basename(colnames(df.i)[ncol(df.i)]),".bam",1)
  colnames(df.i)[ncol(df.i)] <- sample 
  df.i2 <- df.i[c("Geneid",sample)]
  
  # test duplicated gene id
  test <- table(duplicated(df.i2$Geneid))
  if(length(test) != 1){
    LuckyVerbose(sample,": duplicated id existed!")
    colnames(df.i2)[2] <- "merge"
    df.i3 <- ddply(df.i2,"Geneid",summarise,counts=sum(merge))
    colnames(df.i3)[2] <- sample
  } else {
    LuckyVerbose(sample,": no duplicated id.")
    df.i3 <- df.i2
  }
  
  # get unique id
  ID <- unique(c(ID,as.character(df.i3$Geneid)))
  
  # get list
  l[[sample]] <- df.i3
}

## merge counts matrix
df <- data.frame(Geneid = ID)
for(i in 1:length(l)){ # i=1
  df.i <- l[[i]]
  df <- left_join(df,df.i,by="Geneid")
}

## tidy up the matrix 
df <- as.matrix(df)
df[is.na(df)] <- 0 #将未检测到的基因的read数计为0
rownames(df) <- Fastextra(as.character(df[,"Geneid"]),"[.]",1)
COUNT <- df[,-match("Geneid",colnames(df))]

## save data
COUNT <- matrix(as.numeric(COUNT),nrow = nrow(COUNT),dimnames = list(rownames(COUNT),colnames(COUNT)))
saveRDS(COUNT,file=paste0(path_count,"/matrix_counts.rds"))
write.table(COUNT,paste0(path_count,"/matrix_counts.txt"),sep = "\t",col.names = T,row.names = T,quote = F)

LuckyVerbose("All done!")



