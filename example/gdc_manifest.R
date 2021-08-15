

##=========================Description====================##
# fast way to get small RNA counts from SPAR result

##==========================Input========================##
# path_manifest='/data/PanCan/PanCanAtlas-Splicing-2018/manifest.txt'
args <- commandArgs(T)
path_manifest <- args[1]

## Note:
# id	filename	md5	size
# 0de8f3ce-b419-46b7-9983-fa77d2818cd4	merge_graphs_mutex_exons_C2.gff3.gz	c3918d108f5a74d5293c16bbf0ecb787	118710490
# 104c3fbf-ccd7-4d2c-a422-d9d44644fcb0	gtex_merge_graphs_alt_3prime_C2.confirmed.txt.gz	f04e25884645524306e1264fe8b29ddb	201340892

##===========================Usage======================##
# nohup gdc_manifest /data/PanCan/PanCanAtlas-Splicing-2018/manifest.txt &



##==============assistant base function=================##
sc <- "/home/huangwb8/bin/ShellBase.R"
source(sc, encoding = "UTF-8")
LuckyVerbose("Load source functions from: ",sc)

##====================environment=======================##
nd.pac=c("dplyr","plyr","readr")
Plus.library(nd.pac)

##===========================Programe=======================##

a <- read.table(path_manifest,header = T)

id <- as.character(a$id)
name <- as.character(a$filename)
path_ws <- Fastextra(path_manifest,basename(path_manifest),1)

## download with api of GDC
for(i in 1:length(id)){ # i=1
  id.i <- id[i]
  n.i <- name[i]
  url.i <- paste0("https://api.gdc.cancer.gov/data/",id.i)
  e.i <- paste0("wget -bc -T 120 -t inf -P ",path_ws," -o ",n.i,".log -O ",n.i," ",url.i)
  system(e.i)
}
