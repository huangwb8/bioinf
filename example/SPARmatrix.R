

##=========================Description====================##
# fast way to get small RNA counts from SPAR result

##==========================Input========================##
# path_spar='/home/huangwb8/Test/smallRNA-Seq_2/output/SPAR'
args <- commandArgs(T)
path_spar = args[1]
path_save=paste0(path_spar,"/SPARmatrix"); dir.create(path_save,showWarnings = F,recursive = T)

##==============assistant base function=================##
sc <- "/home/huangwb8/bin/ShellBase.R"
source(sc, encoding = "UTF-8")
LuckyVerbose("Load source functions from: ",sc)
LuckyVerbose("Path-SPAR: ",path_spar)

##====================environment=======================##
nd.pac=c("dplyr","plyr","readxl","tidyr")
scipen = 1

LuckyVerbose("Load grobal options: scipen=",scipen,", needed package: ",paste0(nd.pac,collapse = ","))
Plus.library(nd.pac)
options(scipen = scipen)


##===========================Programe=======================##

## get gene expression matrixes
path_res <- list.files(path = path_spar,pattern = "_bam$",full.names = T)


## tidy matrixes
matrix <- NULL; annot <- NULL
for(i in 1:length(path_res)){ # i=1
  
  path_res_i <- path_res[i]
  ge <- list.files(path = path_res_i,pattern = "smRNA_gene_expression.xls", include.dirs = T,recursive = T,full.names = T)

  
  ## create matrix.i for every align type
  matrix.i <- NULL; annot.i <- NULL
  for( j in 1:length(ge)){ # j=1
    ge.j <- ge[j]
    
    res.j <- read.table(ge.j,sep = "\t",comment.char = "",header = T,check.names = F,stringsAsFactors = F)
    type.j <- Fastextra(ge.j,"[/]") %>% .[grep("_bam",.)] %>% Fastextra(.,"_bam",1) # colnames(res.j)
    matrix.j <- res.j[-match("GeneClass",colnames(res.j))]
    base.j <- Fastextra(ge.j,"[/]") %>% .[length(.)-2]
    colnames(matrix.j)[1] <- c("id")
    x <- colnames(matrix.j)[2:ncol(matrix.j)]
    colnames(matrix.j)[2:ncol(matrix.j)] <- paste(type.j,base.j,x,sep = "_")
    if(is.null(matrix.i)){
      matrix.i <- matrix.j
    } else {
      matrix.i <- left_join(matrix.i,matrix.j,by="id")
    }
    
    ## annot
    annot.j <- res.j[c("#Gene","GeneClass")]
    annot.i <- rbind(annot.i,annot.j)
  }
  
  ## merge matrix.i
  if(is.null(matrix)){
    
    matrix <- matrix.i
  } else {
    matrix <- left_join(matrix,matrix.i,by="id")
  }

  ## annot
  annot <- rbind(annot,annot.i)
  
  
}

## remove character Aligned.sortedByCoord.out
colnames(matrix) <- gsub("Aligned.sortedByCoord.out","",colnames(matrix))

## merge duplicates
sm_rna <- Fastextra(as.character(matrix$id),'[+|-][:]',2)

mt <- mergeMatrixDup(x = matrix[2:ncol(matrix)], 
                            mergeCol = F,
                            mergeRow = T, 
                            fun_row = sum, 
                            refRow = sm_rna)

## annotation
LuckyVerbose("Annotation...")
annot[,1] <- Fastextra(as.character(annot[,1]),'[+|-][:]',2)
annot <- annot[!duplicated(as.character(annot[,1])),]
rownames(annot) <- as.character(annot[,1])
annot0 <- annot[rownames(mt),]
colnames(annot0) <- gsub("#","",colnames(annot0))
LuckyVerbose("Annotation done!")

## merge
mt2 <- cbind(annot0,mt)
saveRDS(mt2,file = paste0(path_save,"/SPARMergeMatrix.rds"))
write.csv(mt2,paste0(path_save,"/SPARMergeMatrix.csv"),row.names = F)

## split
align_type <- colnames(mt) %>% Fastextra(.,"_",1) %>% unique
matrix_type <- c("ReadCount","RPM")
for(i in align_type){# i = "STAR"
  for(j in matrix_type){ # j = "RPM"
    nc <- intersect(grep(i,colnames(mt)),grep(j,colnames(mt)))
    mt.x <- cbind(annot0,mt[nc])
    colnames(mt.x) <- gsub(paste0(i,"_"),"",colnames(mt.x))
    colnames(mt.x) <- gsub(paste0("_",j),"",colnames(mt.x))
    saveRDS(mt.x,file = paste0(path_save,"/",i,"_",j,"Data.rds"))
    write.csv(mt.x,paste0(path_save,"/",i,"_",j,"Data.csv"),row.names = F)
  }
  
}

LuckyVerbose("SPARmatrix: Done!")







