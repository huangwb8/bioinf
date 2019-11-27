

## Get downloading links

### input file is like
## @downloading
## html1
## html2 # some annotation
## ...
## (blank line is available. It would be ignored when loading)
## @downloaded
## html1
## html2 # some annotation
## ...


args <- commandArgs(T) #print(args[1])
dlist <- read.table(args[1],sep = "\t",header = F,check.names = T,comment.char = "#")

# dlist <- read.table("/data/annovar/version_20191024/humandb/hg19_download/hg19_annovar",sep = "\t",header = F,check.names = T,comment.char = "#")


dlist <- as.character(dlist[,1])

line_comment <- grep("@",dlist)
line_downloading <- intersect(line_comment,grep("downloading",dlist))
line_downloaded <- intersect(line_comment,grep("downloaded",dlist))


if(line_downloaded >= line_downloading){
  dlinks <- dlist[setdiff(line_downloading:line_downloaded,
                          c(line_downloading,line_downloaded))]
} else {
  dlinks <- dlist[setdiff(line_downloading:length(dlist),
                          c(line_downloading))]
}
dlinks <- gsub(" ","",dlinks)

## output
for( i in dlinks){
  e.i <- paste0("echo '",i,"'")
  system(e.i)
}





