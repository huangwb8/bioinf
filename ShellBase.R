
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
LuckyVerbose <- function(...,levels = 1,type = NULL, show.sys.time = T){
  ## Verbose type
  if(is.null(type)){
    if(levels == 1){
      type <-  "message"
    } else {
      type <-  "cat"
    }
  }

  ## level symbol
  if(levels > 1){
    ls <- paste(rep("o",(levels-1)),collapse = "")
  } else {
    ls <- ""
  }

  ## System time
  if(show.sys.time){
    ls <- paste0(as.character(format(Sys.time(), "%Y-%m-%d %H:%M:%S")),' | ', ls)
  }

  ## do Verbose
  if(type == "message"){
    return(base::message(ls,...))
  } else if(type == "cat") {
    return(base::cat(ls,...,"\n"))
  } else {
    print("Please input right type.")
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

## lucky::mergeMatrixDup
mergeMatrixDup <- function(x,
                           mergeCol = T,
                           fun_col = mean,
                           refCol = NULL,
                           mergeRow = T,
                           fun_row = mean,
                           refRow = NULL,
                           parallel = T){
  
  ## select duplicate data for row
  if(mergeRow & !is.null(refRow)){
    LuckyVerbose("Merge duplicate for row...")
    dupID_row <- refRow[duplicated(refRow)]
    logi <- refRow %in% dupID_row
    x2 <- x[logi, ] # View(x[!logi, ])
    ref <- refRow[grep(T,logi)]
    x2 <- apply(x2,2,function(x)tapply(x,ref,fun_row))
    merN <- rownames(x2)
    x2 <- rbind(x[!logi, ],x2)
    rownames(x2)<- c(refRow[!logi],merN)
  } else {
    LuckyVerbose("Ignore duplicate for row.")
    x2 <- x
  }
  
  ## select duplicate data for col
  if(mergeCol){
    LuckyVerbose("Merge duplicate for col...")
    ## new data
    x <- x2
    
    ## refCol
    if(is.null(refCol)) refCol = colnames(x)
    
    ## select duplicate data
    dupID_Col <- refCol[duplicated(refCol)]
    logi <- refCol %in% dupID_Col
    x2 <- x[,logi] # View(x[,!logi])
    ref <- refCol[grep(T,logi)]
    x2 <- t(apply(x2,1,function(x)tapply(x,ref,fun_col)))
    merN <- colnames(x2)
    x2 <- cbind(x[,!logi],x2)
    colnames(x2) <- c(refCol[!logi],merN)
    
  } else {
    LuckyVerbose("Ignore duplicate for col.")
  }
  
  ## Output data
  LuckyVerbose("All done!")
  return(x2)
  
}




#' @title Is one True?
#' @description Is one True?
#' @param vt a vector
#' @return logic
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' vt <- c(T,"2","1","3");is.one.true(vt)
#' vt <- c(F,"2","1","3");is.one.true(vt)
#' vt <- c(T,2,1,3);is.one.true(vt)
#' vt <- c(T,F);is.one.true(vt)
#' vt <- c(F,F);is.one.true(vt)
#' @export
is.one.true <- function(vt){
 if(!is.logical(vt)){
   #vt不全是逻辑值
   if(!is.numeric(vt)){
     #vt不是数值型变量
     l1 <- T %in% vt
     return(l1)
   } else {
     #vt是数值型向量
     print("向量是数值型，无法区分。")
   }
 } else {
   #vt是逻辑型向量
   l1 <- T %in% vt
   return(l1)
 }

}
