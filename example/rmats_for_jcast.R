

##=========================Input========================##

path_rmats_res = c(
  '~/Project/RNA-Seq_A2B1-KO/output/rmats/STAR_Ensembl_104_unmasked/DM.vs.DB',
  '~/Project/RNA-Seq_A2B1-KO/output/rmats/STAR_Ensembl_104_unmasked/DM.vs.DC'
)
path_jcast_res = c(
  '~/Project/RNA-Seq_A2B1-KO/output/jcast/rmats/STAR_Ensembl_104_unmasked/DM.vs.DB',
  '~/Project/RNA-Seq_A2B1-KO/output/jcast/rmats/STAR_Ensembl_104_unmasked/DM.vs.DC'
)
fdr.cutoff = 0.05


##==============assistant base function=================##
source("/home/huangwb8/bin/ShellBase.R",encoding = "UTF-8")


##====================environment=======================##
nd.pac=c("dplyr","plyr","readr")
Plus.library(nd.pac)

# grobal options:
options(scipen = 100)

##===========================Programe====================##


for(p in 1:length(path_rmats_res)){ # p=1
  
  path_rmats_res.p <- path_rmats_res[p]
  path_jcast_res.p <- path_jcast_res[p]
  
  dir.create(paste0(path_jcast_res.p,'/rmats_res/'),showWarnings = F, recursive = T)
  
  f <- list.files(path = path_rmats_res.p, pattern = 'JCEC.txt$|JC.txt$',full.names = T,recursive = F)
  
  for(i in 1:length(f)){ # i=1
    
    f.i <- f[i]
    n.i <- rev(Fastextra(f.i,'[/]'))[1]
    df.i <- read_delim(f.i, delim = '\t') %>% filter(FDR < fdr.cutoff, !is.na(geneSymbol))
    # df.i <- read.table(f.i, sep = '\t',header = T)
    write.table(df.i, paste0(path_jcast_res.p,'/rmats_res/',n.i), sep = '\t', quote = F,row.names = F)
    
  }
  
}





# Test
if(F){
  
  df.i <- read_delim(f.i, delim = '\t')
  df.i.na <- df.i[is.na(df.i$geneSymbol),]

}








