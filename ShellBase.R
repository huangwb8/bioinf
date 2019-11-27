
## ShellBase.R
# source("/home/huangwb8/bin/ShellBase.R",encoding = "UTF-8")


## Description:
# Basic functions useful and used-frequently


## ?lucky::Fastextra
Fastextra <- function(vt, 
                       split, 
                       n = NULL) {
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

## ?lucky::Fastmatch
Fastmatch <- function(pattern, x){
  p <- apply(as.matrix(pattern), 1, function(z) match(z, x))
  return(p)
}


## ?lucky::Fastgrep
Fastgrep <- function(pattern, x){
  p <- apply(as.matrix(pattern), 1, function(z) grep(z, x))
  p <- unlist(p)
  p <- as.vector(as.matrix(p))
  return(p)
}

## ?lucky::Plus.library
Plus.library <- function(packages){
  x <- installed.packages()
  x <- as.data.frame(x)
  installed <- as.character(x$Package)
  is <- length(grep("pacman", installed))
  if (is == 0) {
    install.packages("pacman")
  }
  
  ## load package
  for (i in packages) {
    pacman::p_load(char = i)
  }
}


## ?lucky::LuckyVerbose
LuckyVerbose <- function(..., 
                         levels = 1, 
                         type = NULL){
  if (is.null(type)) {
    if (levels == 1) {
      type <- "message"
    }
    else {
      type <- "cat"
    }
  }
  if (levels > 1) {
    s1 <- paste(rep(" ", (levels - 1)), collapse = "")
    s2 <- paste(rep("o", (levels - 1)), collapse = "")
    ls <- paste(s1, s2, collapse = "")
  }
  else {
    ls <- ""
  }
  if (type == "message") {
    return(base::message(ls, " ", ...))
  }
  else {
    if (type == "cat") {
      return(base::cat(ls, ..., "\n"))
    }
    else {
      print("Input right type.")
    }
  }
}

## ?lucky::cut_vector
cut_vector <- function(vt,nsplit=100){
  if(nsplit==1){
    v2 <- list(vt);
    names(v2) <- paste0("1-",length(v2))
  } else {
    ##转化成位置向量
    len.vt=1:length(vt)
    
    ##间隔
    len1 <- floor(length(len.vt)/nsplit)
    
    ## low.ci and upper.ci
    low.ci <- 1;
    for(i in 1:(nsplit-1)){
      low.ci[i+1] <- low.ci[i] + len1
    }
    upper.ci <- len1
    for(i in 1:(nsplit-1)){
      upper.ci[i+1] <- upper.ci[i] + len1
    }
    
    ## upper.ci的最后一个值
    upper.ci[length(upper.ci)] <- length(vt)
    
    ## 形成列表
    v2 <- NULL
    for(i in 1:length(low.ci)){
      v2.i <- low.ci[i]:upper.ci[i]
      v2.i <- vt[v2.i]
      v2 <- c(v2,list(v2.i))
      names(v2)[i] <- paste0(low.ci[i],"-",upper.ci[i])
    }
  }
  ## 输出结果
  return(v2)
}





