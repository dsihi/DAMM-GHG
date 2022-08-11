source("proc.formula.R")
lme.inla.formula <- function(formula,data,param=NULL)
{
  ## helper function to formulate INLA formula
  ## formula: mixed effect formula in lme4 syntax
  ## data : data frame to evaluate the data
  ## param: a list of hyper-parameters consistent with each random terms
  ##      for univariate random, default to (0,0.00001) for truncated normal prior
  ##      for multivariate random, default to identity metrix and df equal 2 plus dimension
  ##      multivariate random can't have more than five terms due to INLA constraints
  ## value: a vector of INLA formula terms
  ##        index data frame useful for model fitting
  ## assumption: only single group, no nesting allowed
  ## 
  ## this function is based on the proc. formula function in marked package
  ##
  lst <- proc.form(formula) ## process mixed formula
  
  if(length(lst$re.model)==0){
    return(
      list(formulas=labels(terms(formula)),index=NULL)
    )
  }
  r <- vector("list",length(lst$re.model)) 
  for(i in 1:length(r)){  ## loop through random terms
    elmt <- lst$re.model[[i]]
    dm <- model.matrix(formula(elmt$model),data)
    v0 <- all.vars(formula(elmt$model))
    g0 <- all.vars(formula(elmt$sub))
    
    ## make copies of index
    ## ASSUMING ONE GROUP ONLY, No NESTING
    if(ncol(dm)==1){
      ind1 <- NULL
    } 
    else{
      ng <- length(unique(data[,g0]))
      ind1 <- matrix(as.numeric(data[,g0]),nrow(data),ncol(dm))
      ind1 <- ind1+rep((0:(ncol(dm)-1))*ng,each=nrow(data))
      colnames(ind1) <- paste(g0,1:ncol(dm),sep="")
    }
    
    ## make formula
    out <- rep("",ncol(dm))
    if(is.null(param)){
      ## default hyper-parameter
      tmp <- diag(ncol(dm))
      parami <- c(ncol(dm)+2,1,1,0)
    }
    else{
      parami <- param[[i]]
    }
    if(ncol(dm)==1){
      out[1] <- paste(
        "f(",g0,",",
        'model="iid",',
        'hyper=list(prec=list(prior="logtnormal",param=c(0,0.00001))))',
        sep="")
    }
    if(ncol(dm)>5){
      stop(ncol(dm)," dimension not supported by INLA\n")
    }
    if(ncol(dm)>1){
      out[1] <- paste(
        "f(",colnames(ind1)[1],",",
        'model="iid',ncol(dm),'d",',
        "n=",ncol(dm)*ng,
        ",hyper=list(theta1=list(param=c(",paste(parami,collapse =","),"))))",
        sep="")
      for(j in 2:ncol(dm)){
        out[j] <- paste("f(",colnames(ind1)[j],",",
                        colnames(dm)[j],',copy="',
                        colnames(ind1)[1],'")',sep="")
      }
    }
    r[[i]] <- list(formula=out,index=ind1)
  }
  
  tmp <- do.call(c,lapply(r,function(elmt) elmt$formula))
  idx <- do.call(cbind,lapply(r,function(elmt) elmt$index))
  
  list(formulas=c(lst$fix.model,tmp),index=idx)
}
test <- function()
{
  rm(list=ls())
  source("lme_inla_formula_2.R")
  source("proc.formula.R")
  df <- data.frame(x=rnorm(100),siteID=1:10,month=rep(1:4,25))
  r0 <- lme.inla.formula(~(x|siteID)+(1|month),data=df)
  r0$formulas
}

test2 <- function()
{
  rm(list=ls())
  source("lme_inla_formula_2.R")
  source("proc.formula.R")
  df <- data.frame(x=rnorm(100),siteID=1:10,month=rep(1:4,25))
  r0 <- lme.inla.formula(x~month,data=df)
  r0$formulas
}

develop <- function()
{
  source("proc.formula.R")
  lst <- proc.form(~(x|siteID)+(1|month))
  param <- NULL
  r <- vector("list",length(lst$re.model))
  for(i in 1:length(r)){
    elmt <- lst$re.model[[i]]
    dm <- model.matrix(formula(elmt$model),df)
    v0 <- all.vars(formula(elmt$model))
    g0 <- all.vars(formula(elmt$sub))
    
    ## make copies of index
    ## ASSUMING ONE GROUP ONLY, No NESTING
    ng <- length(unique(df[,g0]))
    if(ncol(dm)==1){
      ind1 <- NULL
    } 
    else{
      ind1 <- matrix(as.numeric(df[,g0]),nrow(df),ncol(dm))
      ind1 <- ind1+rep((0:(ncol(dm)-1))*ng,each=nrow(df))
      colnames(ind1) <- paste(g0,1:ncol(dm),sep="")
    }
    
    ## make formula
    out <- rep("",ncol(dm))
    if(is.null(param)){
      tmp <- diag(ncol(dm))
      param <- c(ncol(dm)+2,diag(tmp),tmp[upper.tri(tmp)])
    }
    if(ncol(dm)==1){
      out[1] <- paste(
        "f(",g0,",",
        'model="iid",',
        'hyper=list(prec=list(prior="logtnormal",param=c(0,0.00001))))',
        sep="")
      
    }
    if(ncol(dm)>5){
      stop(ncol(dm)," dimension not supported by INLA\n")
    }
    if(ncol(dm)>1){
      out[1] <- paste(
        "f(",colnames(ind1)[1],",",
        'model="iid',ncol(dm),'d",',
        "n=",ncol(dm)*ng,
        ",hyper=list(theta1=list(param=c(",paste(param,collapse =","),"))))",
        sep="")
      for(j in 2:ncol(dm)){
        out[j] <- paste("f(",colnames(ind1)[j],",",
                        colnames(dm)[j],',copy="',
                        colnames(ind1)[1],'")',sep="")
      }
    }
    r[[i]] <- list(formula=out,index=ind1)
  }
  
  tmp <- do.call(c,lapply(r,function(elmt) elmt$formula))
  idx <- do.call(cbind,lapply(r,function(elmt) elmt$index))
}

