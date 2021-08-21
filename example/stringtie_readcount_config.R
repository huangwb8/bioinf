

##==========================Input========================##

path_stringtie_res = '~/Project/RNA-Seq_A2B1-KO/output/stringtie/hisat2'



##==============assistant base function=================##
source("/home/huangwb8/bin/ShellBase.R",encoding = "UTF-8")


##====================environment=======================##
nd.pac=c("dplyr","plyr","readr","stringi","tidyr")
Plus.library(nd.pac)

# grobal options:
options(scipen = 100)

##===========================Programe====================##

f <- list.files(path = path_stringtie_res, pattern = 'merged.gtf$', full.names = T, recursive = T) %>% stringi::stri_reverse() 

f.name <- Fastextra(f, '[/]', 2) %>% stringi::stri_reverse() 

df.f <- cbind(name = f.name, path = stri_reverse(f))

write.table(df.f, paste0(path_stringtie_res, '/preDE.config.txt'), sep = '\t', row.names = F,col.names = F, quote = F)






