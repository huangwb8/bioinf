



##==========================Input========================##
# process=32
# path_interval_list="E:/@Analysis/gatk/gatk_internal.list_hg38.exon.bed"
# path_scatter_interval="E:/@Analysis/gatk/scatter_interval"
# bed_type="exon_hg38"
replace = T

## Note:
# The content of interval_list should be like: 
# X1        X2     X3
# <chr>  <int>  <int>
# chr1   69090  70008
# chr1  450739 451678
# chr1  685715 686654
# chr1  925941 926013
# chr1  930154 930336
# chr1  931038 931089

##==============assistant base function=================##
sc <- "/home/huangwb8/bin/ShellBase.R"
source(sc, encoding = "UTF-8")
LuckyVerbose("Load source functions from: ",sc)

##====================environment=======================##
nd.pac=c("readr")
scipen = 1

LuckyVerbose("Load grobal options: scipen=",scipen,", needed package: ",paste0(nd.pac,collapse = ","))
Plus.library(nd.pac)
options(scipen = scipen)


##=========================Programe====================##

test <- list.files(path_scatter_interval)
if(length(test)>0 & !replace){
  paste0("Scatter interval file had been existed:   ",path_scatter_interval)
  list.files(path_scatter_interval)
} else {
  
  ## create new files or replace old files
  ## Get bed data
  bed <- read_delim(path_interval_list,delim = "\t",col_names = F) # str(bed)
  colnames(bed) <- c("chr","start","end")
  bed$interval <- bed$end - bed$start
  sum_interval <- sum(bed$interval)
  avr_interval <- floor(sum_interval/process)
  
  ## Calculate accumulative sum for every loci
  LuckyVerbose("Calculate accumulative sum for every loci...")
  accumulate_sum <- NULL
  for(i in 1:nrow(bed)){
    accumulate_sum[i] <- sum(bed$interval[1:i])
  }
  
  ## Splite bed via accumulative base length
  LuckyVerbose("Get cut off position...")
  cut_off_position <- NULL
  for(i in 1:process){ # i=1
    cut_off_position[i] <- sum(accumulate_sum < i*avr_interval)
  }
  new_cut_off <- cut_off_position[-process]
  new_cut_off[process] <- nrow(bed)
  new_cut_off <- c(0,new_cut_off)
  # View(bed[new_cut_off,])
  
  ## save splited bed
  dir.create(path_scatter_interval,
             recursive = T,
             showWarnings = F)
  file.remove(list.files(path_scatter_interval,full.names = T))
  for(i in 1:(length(new_cut_off)-1)){ # i=1
    bed.i <- bed[(new_cut_off[i]+1):new_cut_off[i+1],]
    order.i <- ifelse(nchar(as.character(i))==1,paste0("00",i),ifelse(nchar(as.character(i))==2,paste0("0",i),i))
    name.i <- paste0(path_scatter_interval,"/",bed_type,"_",order.i,"_",(new_cut_off[i]+1),"-",new_cut_off[i+1],".bed")
    write.table(bed.i[1:3],file = name.i,sep = "\t",quote = F,col.names = F,row.names = F)
  }
  
  ## End function
  list.files(path_scatter_interval)
  LuckyVerbose("All done!")
}


