##############################################################################
# Script information
# Title: Transform the pseudo-bulk expression matrix
# Author: Angli Xue
# Date: 2022-05-17
# Description: This R script was written to transform, scale, standardize pseudo-bulk expression matrix
##############################################################################

pseudobulk_scaling <- function(expr=NA,pi0=F,log1p=F,scale=F,RINT=F,HVG=F){
  if(all( isFALSE(pi0) & isFALSE(log1p) & isFALSE(scale) & isFALSE(RINT) & isFALSE(HVG) )){
    message("No filter or transformation selected! Return the input matrix.")  
    return(expr)
  }
  message("Transforming the pseudo-bulk expression matrix...")  
  # Check data type of the input of expr
  if(is.data.frame(expr) | is.matrix(expr)){
    expr = as.data.frame(expr)
  } else{
    stop("Please input a data frame or matrix.")
  }
  
  message(paste0("The input matrix has ",nrow(expr)," rows and ",ncol(expr)," cols"))  
  # Check NAs in the matrix
  idx=which(is.na(rowMeans(expr)))
  if(length(idx) > 0){
  warning(paste0(length(idx)," rows have NA values"))
  expr = expr[-idx,]
  }
  if(nrow(expr) > ncol(expr)){
  warning("The number of genes is smaller than the sample size. Please double check your data!")
  }
  message(paste0("The settings are:\n","pi0 = ",pi0,"\nlog1p = ",log1p), "\nscale = ",scale, "\nRINT = ", RINT, "\nHVG = ", HVG)
  
  ## Exclude zero-inflated genes
  if(isFALSE(pi0)){
    message("No filter for pi0")
    max_var = ncol(expr)
    expr2 = expr
  } else if(is.numeric(pi0) & pi0>=0 & pi0<=1){
    # Calculate the pi0 for each variable
    pi_hat=colSums(expr==0,na.rm=T)/sum(!is.na(expr[,1]))
    if(sum(pi_hat >= pi0)==0){
      message("No genes with pi0 >= ",pi0)
      expr2 = expr
    } else{
      message(paste0("Excluded ",sum(pi_hat >= pi0)," genes with pi0 >= ",pi0))
      expr2=expr[,which(pi_hat < pi0)]
    }
  } else{
    stop("pi0 out of bound or wrong data type! Please specificy a number between 0 and 1.")
    }
  
  ## Whether use top highly variable genes (HVGs)
  if(!exists("max_var")){max_var = ncol(expr2)}
  if(isFALSE(HVG)){
    message(paste0("Use all ",max_var," remaining genes"))
  } else if(is.numeric(HVG) & HVG > 0 & HVG <= max_var){
    message("Use top ",HVG," HVGs to generate PEER factors")

    # Avoid genes with extremely low expression (i.e., < 0.001)
    cm=as.numeric(colMeans(expr2,na.rm=T))
    if(sum(cm >= 0.001) > 0){
      expr2=expr2[,which(cm>=0.001)]
    }
    if(ncol(expr2) < HVG){
      stop(paste0("Too many genes with extremely low expression. Please set the number of HVG lower than ",ncol(expr2),"."))
    }
 
    # Calculate the mean, variance, and CV
    cm=as.numeric(colMeans(expr2,na.rm=T))
    cv=apply(expr2,2,function(x)var(x,na.rm=T))/cm
    cv=as.numeric(cv)
    expr2=expr2[,which(cv>=sort(cv,decreasing=T)[HVG])]
  } else{
    stop("HVG out of bound or wrong data type! Please specificy a number between 0 and total number of variables (",max_var,").")
  }

  ## Transformation and standardization
  if((!RINT & !log1p & !scale)){
    message("No transformation or scaling is selected")
    expr2=as.data.frame(expr2) 
    message("Transformation or scaling completed!")  
    return(expr2)
  }
  
  if(RINT){
    require(RNOmni)
    message("Use Rank-based Inverse Normal Transformation (RINT)")
    message("The log1p and scale flags will be turned off")
    # RINT
    expr2=apply(expr2,2,function(x)RankNorm(x))
    # Scale and center
    expr2=scale(expr2, scale = T,center = T)
    expr2=as.data.frame(expr2) 
    message("Transformation or scaling completed!")  
    return(expr2)
  } 
  
  if(log1p){
    if(sum(expr2<0) > 0){
    stop("Negative values found in matrix. Can not use log transform. Please check your data!")
    }
    message("Use log(x+1) transformation")
    # Log(x+1) transformation
    expr2=apply(expr2,2,function(x)log(x+1))
  } 
  
  if(scale){
    message("Standardize the matrix column-wise (scale and center)")
    # Convert each column to mean = 0 and sd = 1
    expr2=scale(expr2, scale = T,center = T)
  } 
  
  # Make sure the output is in data.frame
  expr2=as.data.frame(expr2) 
  message("Transformation and scaling completed!")  
  return(expr2)
}

####
