
##=========================Input========================##

path_stringtie_res = '~/Project/RNA-Seq_A2B1-KO/output/stringtie/hisat2'
path_suppa_res = '~/Project/RNA-Seq_A2B1-KO/output/suppa/stringtie/hisat2'


# Sample & Group
Group = list(
  'DM' = c('D1','D3','D8'),
  'DB' = c('D4','D5','D9'),
  'DC' = c('D6','D2','D7')
)

event.type = c('A3','A5','AF','AL', 'MX', 'RI', 'SE')

genome_type = 'hg38'

##==============assistant base function=================##
source("/home/huangwb8/bin/ShellBase.R",encoding = "UTF-8")

##====================environment=======================##
nd.pac=c("dplyr","plyr","readr")
Plus.library(nd.pac)

# grobal options:
options(scipen = 100)

##===========================Programe====================##

# TPM matrix
if(T){
  
  f <- list.files(path = path_stringtie_res, pattern = '.transcript.tpm.tab',full.names = T,recursive = T)
  for(i in 1:length(Group)){ # i=1
    
    g.i <- Group[[i]]
    n.i <- names(Group)[i]
    f.i <- f[grepl(paste0(g.i,collapse = '|'),f)]
    df.i <- NULL
    for(j in 1:length(f.i)){ # j=1
      
      df.j <- read.table(f.i[j])
      df.j$transcript <- rownames(df.j)
      if(is.null(df.i)){
        df.i <- df.j
      } else {
        df.i <- full_join(df.i,df.j,by='transcript')
      }
      
    }
    rownames(df.i) <- as.character(df.i$transcript)
    df.i <- df.i[-match('transcript',colnames(df.i))]
    for(j in 1:ncol(df.i)) df.i[,j][is.na(df.i[,j])] <- 0
    
    write.table(df.i, paste0(path_suppa_res,'/', n.i, '.stringtie.tpm'),sep = '\t',quote = F)
  }
  
}

# PSI: Differential transcript usage
if(T){
  
  f <- list.files(path = path_suppa_res, pattern = paste0('_isoform.psi'), full.names = T,recursive = T) %>% .[grepl(., pattern = genome_type)]
  psi_transcript <- list()
  for(g in 1:length(Group)){ # g=1
    name.g <- names(Group)[g]
    sample.g <- Group[[g]]
    f.g <- f[grepl(paste0(sample.g , collapse = '|'),f)]
    df.g <- NULL
    for(j in 1:length(f.g)){ # j=1
      f.j <- f.g[j]
      df.j <- read.table(f.j)
      df.j$event <- rownames(df.j)
      if(is.null(df.g)){
        df.g <- df.j
      } else {
        df.g <- full_join(df.g,df.j,by = 'event')
      }
    }
    rownames(df.g) <- as.character(df.g$event)
    df.g <- df.g[-match('event',colnames(df.g))]
    for(j in 1:ncol(df.g)){ # j=1
      df.g[,j] <- as.character(df.g[,j])
      df.g[,j][df.g[,j] %in% c('NaN',NaN)] <- 'nan'
      df.g[,j][df.g[,j] %in% c('NA',NA)] <- 'nan'
      df.g[,j][is.na(df.g[,j])] <- 'nan'
      # df.g[,j] <- as.numeric(df.g[,j])
      # df.g[,j][df.g[,j] %in% c('NaN',NaN)] <- 0
      # df.g[,j][df.g[,j] %in% c('NA',NA)] <- 0
      # df.g[,j][is.na(df.g[,j])] <- 0
    }
    # write.table(df.g,paste0(path_suppa_res,'/', name.g, '_isoform.psi'), sep = '\t', quote = F)
    psi_transcript[[name.g]] <- df.g
  }
  
  # Force to same events
  same_event <- NULL
  for(i in 1:length(psi_transcript)){ # i=1
    df.i <- psi_transcript[[i]]
    event.i <- rownames(df.i)
    if(is.null(same_event)){
      same_event <- event.i
    } else {
      same_event <- intersect(same_event, event.i)
    }
  }
  
  # Save files
  for(g in 1:length(psi_transcript)){ # g=1
    df.g <- psi_transcript[[g]][same_event,]
    name.g <- names(psi_transcript)[g]
    write.table(df.g, paste0(path_suppa_res,'/', name.g, '_isoform.psi'), sep = '\t', quote = F)
  }

}

# PSI: Differential splicing with local events
if(T){
  
  psi_local <- list()
  for(i in 1:length(event.type)){ # i=1
    
    et.i <- event.type[i]
    f.i <- list.files(path = path_suppa_res, pattern = paste0('_',et.i,'.psi'),full.names = T,recursive = T) %>% .[grepl(., pattern = genome_type)]
    
    for(g in 1:length(Group)){ # g=1
      
      name.g <- names(Group)[g]
      sample.g <- Group[[g]]
      f.g <- f.i[grepl(paste0(sample.g , collapse = '|'),f.i)]
      
      df.g <- NULL
      for(j in 1:length(f.g)){ # j=1
        
        f.j <- f.g[j]
        df.j <- read.table(f.j)
        df.j$event <- rownames(df.j)
        if(is.null(df.g)){
          df.g <- df.j
        } else {
          df.g <- full_join(df.g,df.j,by = 'event')
        }
      }
      rownames(df.g) <- as.character(df.g$event)
      df.g <- df.g[-match('event',colnames(df.g))]
      for(j in 1:ncol(df.g)){ # j=1
        df.g[,j] <- as.character(df.g[,j])
        df.g[,j][df.g[,j] %in% c('NaN',NaN)] <- 'nan'
        df.g[,j][df.g[,j] %in% c('NA',NA)] <- 'nan'
        df.g[,j][is.na(df.g[,j])] <- 'nan'
        # df.g[,j] <- as.numeric(df.g[,j])
        # df.g[,j][df.g[,j] %in% c('NaN',NaN)] <- 0
        # df.g[,j][df.g[,j] %in% c('NA',NA)] <- 0
        # df.g[,j][is.na(df.g[,j])] <- 0
      }
      psi_local[[et.i]][[name.g]] <- df.g
      # write.table(df.g,paste0(path_suppa_res,'/',name.g,'.',et.i,'.psi'),sep = '\t',quote = F)
      
    }
    
  }
  
  # Force to same events
  same_event <- NULL
  for(i in 1:length(psi_local)){ # i=1
    
    event.i <- names(psi_local)[i]
    psi_local.i <- psi_local[[i]]
    
    same_event.i <- NULL
    for(j in 1:length(psi_local.i)){ # j=1
      df.j <- psi_local.i[[j]]
      same_event.j <- rownames(df.j)
      if(is.null(same_event.i)){
        same_event.i <- same_event.j
      } else {
        same_event.i <- intersect(same_event.i, same_event.j)
      }
    }
    
    # Save files
    for(g in 1:length(psi_local.i)){ # g=1
      df.g <- psi_local.i[[g]][same_event.i,]
      name.g <- names(psi_local.i)[g]
      write.table(df.g, paste0(path_suppa_res,'/', name.g, '.', event.i, '.psi'), sep = '\t', quote = F)
    }

  }

}



