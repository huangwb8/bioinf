

## Description
# Function like lucky::Fastextra in R


## example
# Rscript fast_extra.R 'annovar refGene;cytoBand;annovar exac03;annovar avsnp147;annovar dbnsfp30a' '[;]'
# echo "'$database'" | xargs -i fast_extra {} '[;]'

## Programe
args <- commandArgs(T)


Fastextra <- function(vt,split,n = NULL) {
  vt <- as.character(vt)
  get1 <- function(i, split, n = NULL) {
    if (is.null(n)) {
      vt1.i <- unlist(strsplit(i, split))
    }
    else {
      vt1.i <- unlist(strsplit(i, split))[n]
    }
    return(vt1.i)
  }
  vt1 <- apply(as.matrix(vt), 1, function(z) get1(z, split = split, 
                                                  n = n))
  vt1 <- as.vector(vt1)
  return(vt1)
}

## split
a <- Fastextra(args[1],args[2])
for(i in a){
  e.i <- paste0("echo '",i,"'")
  system(e.i)
}



