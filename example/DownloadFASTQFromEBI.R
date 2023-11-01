

#============Information
# Version: 0.1.10
# Author: Weibin Huang
# Download fastq.gz (if any) from https://www.ebi.ac.uk/ena after giving an accession ID like PRJEB402. Support single/pair data.

#============Usage
# In Shell:

# conda activate r4.0

# ws=/home/huangwb8/Test/GSE167977; nohup Rscript /home/huangwb8/Test/GSE167977/bin/DownloadFASTQFromEBI.R --accession PRJNA705699 --path_project $ws --method aspera > $ws/output/fastq/DownloadFASTQFromEBI.log 2>&1 &; tail -f $ws/output/fastq/*.log

# Rscript /home/huangwb8/Test/GSE167977/bin/DownloadFASTQFromEBI.R -h


#============Assistant functions

source("/home/huangwb8/bin/ShellBase.R",encoding = "UTF-8")

# Package
# if(F){
#   # Install "devtools" package
#   if (!requireNamespace("devtools", quietly = TRUE))
#     install.packages("devtools")
#   
#   # Install dependencies
#   if (!requireNamespace("luckyBase", quietly = TRUE))
#     devtools::install_github("huangwb8/luckyBase")
#   
#   # Load packages
#   library(luckyBase)
#   np <- c('curl','optparse'); Plus.library(np)
#   
# }

# Install dependencies
# if (!requireNamespace("curl", quietly = TRUE))
#   install.packages("curl")
# library(curl); 

if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse")
library(optparse);


# Function: downloader
# wget, aspera
fastq_downloader <- function(path_fastq, url, method, losse_check){
  
  # url = "fasp.sra.ebi.ac.uk:/vol1/fastq/SRR138/028/SRR13816328/SRR13816328.fastq.gz"
  url_base = basename(url)

  # Download fastq.gz
  repeat {
    
    # Shell expression
    if(method == 'aspera'){
      # ~/.aspera/connect/bin/ascp -l 1000M -P 33001 -QT -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR138/028/SRR13816328/SRR13816328.fastq.gz ./test.fq.gz
      url_e = paste0('~/.aspera/connect/bin/ascp -l 1000M -P 33001 -QT -k 2 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@', url, ' ', path_fastq,'/',url_base)
      if(losse_check & file.exists(paste0(path_fastq,'/',url_base))){
        url_e = 'completed'
      }
    } else if(method == 'wget'){
      url_e = paste0('wget -c -o ',path_fastq,'/wget_', url_base, '.log -T 120 -t inf --directory-prefix=',path_fastq,' ', url)
    } else {
      url_e <- 'wrongmethod'
    }
    
    
    # Download
    if(url_e == 'completed'){
      
      LuckyVerbose(accession, ': ', url_base, ' had been download completely. Ignored!')
      break
      
    } else if(url_e == 'wrongmethod'){
      
      LuckyVerbose(accession, ": Please use right method. One of wget and aspera.")
      break
      
    } else {
      
      # Download the url
      LuckyVerbose(accession, ': Downloading ', url,'...')
      result <- tryCatch(system(url_e, intern = TRUE), error = function(e) e) # wait = TRUE is default
      error_msg <- as.character(unlist(result))
      is_repeat <- is.one.true(grepl("Error|Partial|Fail", error_msg, ignore.case = T))
      if(is_repeat){
        LuckyVerbose(accession, ": ascp encountered an error. Retrying...")
        Sys.sleep(120) # Another repeat after 120s
      } else {
        LuckyVerbose(accession, ": ",url_base,"'s ascp result: ")
        print(result)
        LuckyVerbose(accession, ": ",url_base,' completed!')
        break
      }
      
    }
    
  }
  
}


#============Parameters

# Shell to R
if(T){
  
  # Define Shell parameters
  option_list <- list(
    make_option(c("-a", "--accession"), type = "character", default = NULL, help = "The accession ID of a GEO project. "),
    make_option(c("-p", "--path_project"), type = "character", default = NULL, help = "The path of project"),
    make_option(c("-m", "--method"), type = "character", default = 'aspera', help = "One of 'aspera', 'wget'"),
    make_option(c("-l", "--losse_check"), action = "store_true", default = TRUE, help = "Wheter to ignore files with suffix 'fastq.gz'. Only works in aspera method!")
  )
  
  # Shell to R
  opt <- parse_args(OptionParser(option_list = option_list))
  
  # R parameters
  accession = opt$accession
  path_project = opt$path_project
  path_fastq = paste0(path_project,'/', 'output/fastq')
  method = opt$method
  losse_check = opt$losse_check
  
} else {
  
  # R parameters
  accession = 'PRJNA705699'
  path_project = 'E:/PythonCloud/@Ubuntu/test/GSE167977'
  path_fastq = paste0(path_project,'/output/fastq')
  method = 'aspera'
  losse_check = T
}

#============Programe
  
# dir
dir.create(path_project, showWarnings = F, recursive = T)
dir.create(path_fastq, showWarnings = F, recursive = T)

# Download filereport file
file_report <- paste0(path_fastq, '/', 'filereport_read_run_',accession, '.tsv')
# if(!file.exists(file_report)){
#   LuckyVerbose(accession, ': Downloading filereport file completed!')
# } else {
#   LuckyVerbose(accession, ': Filereport file exist. Ignore!')
# }
LuckyVerbose(accession, ': Downloading filereport file...')
# curl_download(paste0('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=',accession,'&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,fastq_aspera,submitted_ftp,sra_ftp,bam_ftp&format=tsv&download=true&limit=0'), paste0(path_project, '/', 'filereport_read_run_',accession, '.tsv'), quiet = TRUE,  mode = "wb")
curl_e <- paste0('curl -C - -o ', file_report,' ', "'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",accession,"&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,fastq_aspera,submitted_ftp,sra_ftp,bam_ftp&format=tsv&download=true&limit=0'")
system(curl_e)

# Get URLs
df <- read.table(file_report, sep = '\t', header = T, check.names = F)

# Download
LuckyVerbose(accession, ': Use ',method,' method! ')
fq_urls <- as.character(df$fastq_aspera)
for(urlP in fq_urls){
  # urlP = "fasp.sra.ebi.ac.uk:/vol1/fastq/SRR588/002/SRR5885322/SRR5885322_1.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/SRR588/002/SRR5885322/SRR5885322_2.fastq.gz"
  urls <- Fastextra(urlP, ';')
  for(url in urls) fastq_downloader(path_fastq, url, method, losse_check)
}

# Game over
LuckyVerbose(accession, ': All done!')

