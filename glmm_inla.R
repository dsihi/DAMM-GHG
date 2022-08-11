source("lme_inla_formula.R")
source("proc.formula.R")
source("getVars.R")
source("s2inla.R")
## ver 4) experimental - additive support

## helper script to fit a generalized linear model using INLA
## arguments: 
##  formula:  lme type formula object 
##  smooth :  smooth terms in mgcv format
##  data   :  data frame to evaluate formula
##  param  : a list of parameters of Wishart prior for each random term
##  family : family for response variable
##  E      : offset INLA term
## note: no nested random effect now
## value: INLA fitted object

glmm.inla <- function(formula,data,newdata=NULL,family="gaussian",E=NULL,Ntrials=NULL,param=NULL,weights=NULL,
                      verbose=TRUE)
{
  #browser()
  ## get response (if any) and predictor variables
  l.vars <- getVars(formula)
  stopifnot(nchar(l.vars)[1]>0)
  if(!is.null(newdata)){
    if(!(l.vars[1] %in% colnames(newdata))){
      newdata <- cbind(newdata,NA)
      colnames(newdata)[ncol(newdata)] <- l.vars[1]
    }else{
      newdata[,l.vars[1]] <- NA
    }
    data <- rbind(data[,l.vars,drop=FALSE],
                  newdata[,l.vars,drop=FALSE])
  }
  r0 <- lme.inla.formula(formula=formula,data=data,param=param)
  if(length(r0$formulas)==0){
    for0 <- reformulate("1",response=all.vars(formula)[1])
  }
  else{
    for0 <- reformulate(r0$formulas,response=all.vars(formula)[1])
  }
  df0 <- data
  if(!is.null(r0$index)){
    df0 <- cbind(data,r0$index)
  }
  
  #browser()
  sout <- s2inla(for0,data=df0)
  assign(sout$zmat,sout$X)
  #term1 <- paste(as.character(for0)[2],as.character(sout$formula)[2],sep="+")
  #for1 <- as.character(for0)
  #for1[2] <- term1
  #for1 <- formula(for1)
  
  result0 <- inla(sout$formula,data=sout$data,
                  family=family,E=E,Ntrials = Ntrials,weights=weights,
                  control.predictor = list(link=1),
                  control.compute=list(dic=TRUE,config=TRUE),verbose=verbose)
  result0$formula <- formula
  result0$DIC <- with(result0$dic,c(Dbar=mean.deviance,pD=p.eff,DIC=dic))
  result0
}

test <- function()
{
  rm(list=ls())
  require(INLA)
  ## Session > Set Working Dir > To Source File Loc
  source("glmm_inla_4.R")
  source("lme_inla_spde_3.R")
  df <- data.frame(x1=rnorm(100),x2=rnorm(100),siteID=1:10,month=rep(1:4,25),y=rnorm(100),d1=rnorm(100))
  df$z <- exp(df$y)
  newdf <- data.frame(x1=rnorm(10),x2=rnorm(10),siteID=5,month=3,d2=rnorm(100))
  r1 <- glmm.inla(z~s(x2)+(x1|siteID)+(1|month),data=df,newdata=newdf,family="poisson")
  DIC(r1)
}

