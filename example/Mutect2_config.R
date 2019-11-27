

## environment
nd.pac=c("dplyr","tidyr")
x <- as.data.frame(installed.packages())
installed <- as.character(x$Package)
for(i in nd.pac){
  if(! i %in% installed) install.packages(i)
}

## needed input:
# path_sra2fq='/home/huangwb8/Test/wes_2/input/case/sra2case.txt'
# control_group='Norm'

## Note:
## 1. the sra2fq.txt should be like:
# id          name          group sample
# SRR6269853	DCIS_P1_Exon	DCIS	P1
# SRR6269859	DCIS_P10_Exon	DCIS	P10

## Programe:

# read data
x <- read.table(path_sra2fq,stringsAsFactors = F)
colnames(x) <- c("id","name","group","sample")
x2 <- x[c("name","sample","group")]
x3 <- tidyr::spread(x2,group,name)
treat.name <- colnames(x3)[-Fastmatch(c("sample",control_group),colnames(x3))]

# merge 
df <- NULL
for(i in 1:length(treat.name)){ # i=1
  t.i <- treat.name[i]
  df.i <- x3[c("sample",control_group,t.i)]
  colnames(df.i)[match(t.i,colnames(df.i))] <- "treat"
  df <- rbind(df,df.i)
}

df <- dplyr::arrange(df,sample)

## save
path_config <- paste0(gsub(basename(path_sra2fq),"",path_sra2fq),"Mutect2_config")
write.table(df,path_config,sep = " ",col.names = F,row.names = F,quote = F)

## End Mutect2_config.R


