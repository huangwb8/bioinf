



##==========================Input========================##

# path_bam = '~/Project/RNA-Seq_A2B1-KO/output/align/hisat2'
# path_rmats_res = '~/Project/RNA-Seq_A2B1-KO/output/rmats/stringtie/hisat2'
# 
# # Sample & Group
# Group = list(
#   'DM' = c('D1','D3','D8'),
#   'DB' = c('D4','D5','D9'),
#   'DC' = c('D6','D2','D7')
# )


##==============assistant base function=================##
source("/home/huangwb8/bin/ShellBase.R",encoding = "UTF-8")

##====================environment=======================##
nd.pac=c("dplyr","plyr","readr")
Plus.library(nd.pac)

# grobal options:
options(scipen = 100)

##===========================Programe====================##

f <- list.files(path = path_bam, pattern = '.bam$', full.names = T, recursive = F)

for(i in 1:length(Group)){ # i=1
  
  g.i <- Group[[i]]
  n.i <- names(Group)[i]
  path.i <- f[grepl(paste0(g.i,collapse = '|'),f)]
  path.i2 <- paste0(path.i,collapse = ',')
  write.table(path.i2, file = paste0(path_rmats_res,'/',n.i,'.txt'), sep = '\t',quote = F,row.names = F,col.names = F)
  
}


