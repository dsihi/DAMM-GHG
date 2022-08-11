## simulate.R

#install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable", dep=TRUE)


require(INLA)
source("glmm_inla.R")
simulate <- function(formula,data,nsweep=99)
{
  nls.inla <- glmm.inla(formula,data=data,verbose = FALSE)
  tmp <- inla.posterior.sample(n=nsweep,nls.inla)
  samples <- t(sapply(tmp,function(elmt){
    resid <- rnorm(n=nrow(data),mean=0,sd=1/sqrt(elmt$hyperpar[1]))
    elmt$latent[1:nrow(data),1]+resid
  }))
  t(samples)
}

test <- function()
{
  rm(list=ls())
  source("simulate.R")
  trench <- read.csv("2015data.csv",head=T)
  trench$t__ <- 1:nrow(trench)
  Rh.samples <- simulate(Rh~s(t__),data=trench)
  plot(ts(Rh.samples),plot.type="single")
  points(trench$Rh,pch=21,bg="steelblue1",cex=2)
  
  CH4.samples <- simulate(CH4~s(t__),data=trench)
  plot(ts(CH4.samples),plot.type="single")
  points(trench$CH4,pch=21,bg="salmon",cex=2)
  
  apply(CH4.samples,2,mean)
  
  N2O.samples <- simulate(N2O~s(t__),data=trench,nsweep=3)
  plot(ts(N2O.samples),plot.type="single")
  points(trench$N2O,pch=21,bg="salmon",cex=2)
}
