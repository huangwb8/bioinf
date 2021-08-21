


##==========================Input========================##

# path_vast_res = '~/Project/RNA-Seq_A2B1-KO/output/vast-tools'
# 
# # Sample & Group
# Group = list(
#   'DM' = c('D1','D3','D8'),
#   'DB' = c('D4','D5','D9'),
#   'DC' = c('D6','D2','D7')
# )

##==============assistant base function=================##
# source("/home/huangwb8/bin/ShellBase.R",encoding = "UTF-8")


##====================environment=======================##
# nd.pac=c("dplyr","plyr","readr")
# Plus.library(nd.pac)
# 
# # grobal options:
# options(scipen = 100)

##===========================Programe====================##

df <- NULL
for(i in 1:length(Group)){ # i=1
   
  df.i <- data.frame(
    sample = Group[[i]],
    group = names(Group)[i],
    stringsAsFactors = F
  )
  df <- rbind(df,df.i)
  
}

write.table(df, paste0(path_vast_res, '/config_file'), sep = '\t', quote = F, col.names = F,row.names = F)







