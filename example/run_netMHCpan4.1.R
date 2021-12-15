

##=========================Description====================##
# Run netMHCpan in R


##==========================Input========================##

# FASTA file
path_fasta = '/home/huangwb8/Project/RNA-Seq_A2B1-KO/output/jcast/rmats/STAR_Ensembl_104_unmasked/DM.vs.DB/jcast_20210823160407'
pattern_target_fasta = NULL
pattern_excluded_fasta = c('_orphan', 'canonical')

# netMHCpan
path_res = paste0(path_fasta,'/run_netMHCpan4.1'); dir.create(path_res, showWarnings = F, recursive = T)
path_software='/home/huangwb8/Downloads/netMHCpan/4.1'
MHC_species = 'HLA'

# HLA typing
path_hla_type = '~/Project/RNA-Seq_A2B1-KO/input/HLA_TCGA-COAD.rds' # From TSNAdb. Please see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6203688/
# A*02:01  C*12:03 A*03:01 A*24:02 A*01:01 A*11:01 B*07:02 
# 12780    9209    7842    7272    5591    4787    4022

# main parameters
# para_l = '8,9,10,11'


##==============assistant base function=================##
sc <- "/home/huangwb8/bin/ShellBase.R"
source(sc, encoding = "UTF-8")
LuckyVerbose("Load source functions from: ",sc)

##====================environment=======================##
nd.pac=c("dplyr","plyr","readr", 'tidyr')
scipen = 100

LuckyVerbose("Load grobal options: scipen=",scipen,", needed package: ",paste0(nd.pac,collapse = ","))
Plus.library(nd.pac)
options(scipen = scipen)

##===========================Programe===================##

# Fasta files
f <- list.files(path = path_fasta, pattern = '.fasta$', full.names = T, recursive = F)
if(!is.null(pattern_target_fasta)){
  f <- f[grepl(paste0(pattern_target_fasta, collapse = '|'), f)]
}
if(!is.null(pattern_excluded_fasta)){
  f <- f[!grepl(paste0(pattern_excluded_fasta, collapse = '|'), f)]
}

# allele names
df.allele <- read.table(paste0(path_software,'/data/allelenames'), sep = ' ', stringsAsFactors = F)
# Fastextra(df.allele[,1], '-', 1) %>% table()
# BoLA   DLA  Eqca  Gogo  H   HLA    Mamu  Patr   SLA 
# 182    2    1     1     11  11858  415   105    76
df.mhc <- df.allele[grepl(paste0(MHC_species,'-'), as.character(df.allele[,1])),]
mhc.selected <- readRDS(path_hla_type); names(mhc.selected) <- gsub(':','',names(mhc.selected))
df.mhc.selected <- df.mhc[Fastgrep(names(mhc.selected), df.mhc[,2]),] # table(duplicated(df.mhc.selected$V2))
mhc.selected <- unique(df.mhc.selected[,1])
mhc.selected.list <- cut_vector(mhc.selected, nsplit = round(length(mhc.selected)/20))

# Run software
for(i in 1:length(f)){ # i=1
  
  f.i <- f[i]
  n.i <- basename(f.i) %>% Fastextra('.fasta', 1)
  
  for(j in 1:length(mhc.selected.list)){ # j=1
    
    para_a <- paste0(mhc.selected.list[[j]], collapse = ',')
    n.j <- names(mhc.selected.list)[j]
    
    e.j <- paste0(
      paste0(path_software, '/netMHCpan'), ' ',
      paste0('-a ', para_a), ' ',
      paste0(f.i), ' ',
      paste0('> ', path_res, '/', n.i, '.netMHCpan.HLA',n.j,'.out &')
    )
    
    # e.j <- paste0(
    #   paste0(path_software, '/netMHCpan'), ' ',
    #   paste0('-xls'), ' ',
    #   paste0('-a ', para_a), ' ',
    #   paste0(f.i), ' ',
    #   paste0('-xlsfile ', path_res, '/', n.i, '.netMHCpan.HLA',n.j,'.xls &')
    # )
    
    system(e.j)
    
  }
  
}

