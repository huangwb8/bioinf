




##==========================Input========================##
# path_seg='/home/huangwb8/Test/wes_2/output/gatk_cnv_unpair/PoT/CallCopyRatioSegments/DCIS_P10_Exon.called.seg'
args <- commandArgs(T); path_seg <- args[1]


## Note:
# The content of *.called.seg is like:
#CONTIG  START     END     NUM_POINTS_COPY_RATIO  MEAN_LOG2_COPY_RATIO  CALL
#chr10   46807     22551022    1139                 1.097063             +
#chr10   22567601  43646383    754                 -0.116042             0
#chr10   46385593  46391180    2                    -3.521996            -


# The output *.avinput is like:
#Chr	Start	     End      	Ref 	Alt
#13   20797176   21105944   0      - 
#13   20797176   21105944   0      +


##==============assistant base function=================##
sc <- "/home/huangwb8/bin/ShellBase.R"
source(sc, encoding = "UTF-8")
LuckyVerbose("Load source functions from: ",sc)
LuckyVerbose("Load calledSeg2avinput.R")
getCS_one <- function(x){
  
  
  if(F){
    x <- dat.seg[1,]
  }
  
  chr <- Fastextra(as.character(x[,"CONTIG"]),"chr",2)
  start <- as.numeric(x[,"START"])
  end <- as.numeric(x[,"END"])
  ref <- 0
  alt <- 0

  ## data.frame
  d <- data.frame(
    chr=chr,
    start=start,
    end=end,
    ref=ref,
    alt=alt,
    stringsAsFactors = F
  )
  
  ## output
  return(d)
  
  
}

       
##====================environment=======================##
nd.pac=c("dplyr","plyr","readr")
scipen = 1

LuckyVerbose("Load grobal options: scipen=",scipen,", needed package: ",paste0(nd.pac,collapse = ","))
Plus.library(nd.pac)
options(scipen = scipen)

##===========================Programe=======================##
LuckyVerbose("Origin ",path_seg)
dat.seg <- read.table(path_seg,sep = "\t",comment.char = "@",check.names = F,header = T) %>% filter(.,!CALL %in% 0) # colnames(dat.seg)

dat.seg2 <- adply(dat.seg,1,getCS_one) # colnames(dat.seg2)

path_avinput <- Fastextra(path_seg,basename(path_seg),1)
name_avinput <- gsub(".called.seg","",basename(path_seg))

## amplification
dat.seg3 <-  filter(dat.seg2,CALL %in% "+") %>% .[c("chr","start","end","ref","alt")]
write.table(dat.seg3,file = paste0(path_avinput,name_avinput,"_amplification.avinput"),sep = "\t",col.names = F,quote = F,row.names = F)

## deletion
dat.seg3 <-  filter(dat.seg2,CALL %in% "-") %>% .[c("chr","start","end","ref","alt")]
write.table(dat.seg3,file = paste0(path_avinput,name_avinput,"_deletion.avinput"),sep = "\t",col.names = F,quote = F,row.names = F)

LuckyVerbose("Done ",name_avinput,"!")

