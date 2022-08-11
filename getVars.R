## get response variable (if any) and predictor variables from a formula object
getVars <- function(formula)
{
  if(length(formula)<3){
    ## no response
    r <- c("",all.vars(formula))
  }
  else{
    r <- all.vars(formula)
  }
  r
}

test <- function()
{
  rm(list=ls())
  source("getVars_1.R")
  getVars(z~s(d)+(x|siteID)+(1|month))
  getVars(~s(d)+(x|siteID)+(1|month))
}

class(~s(d)+(x|siteID)+(1|month))
length(z~s(d)+(x|siteID)+(1|month))
#length
#one-sided formula will have a length of 2, 
#while the two-sided formula will have a length of 3.