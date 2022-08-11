## damm_utils.R

## utility function for damm

## original scale -> unrestricted
## x : variate in uniform scale
## z : transformed scale on real line
## l : lower limit of the uniform interval
## r : upper limit of the uniform interval
## log: whether to apply log to x 
logit_tran <- function(x,l,r,log=FALSE)
{
  if(log){
    x <- log(x)
    l <- log(l)
    r <- log(r)
  }
  log((x-l)/(r-x))
}
## unrestricted scale -> original
expit_tran <- function(z,l,r,log=FALSE)
{
  if(log){
    l <- log(l)
    r <- log(r)
  }
  l.z <- exp(z)
  output <- ifelse(l.z>0 & is.infinite(l.z),r,(l+r*l.z)/(1+l.z))
  if(log){
    output <- exp(output)
  }
  output
}

## helper function to transform real parameters back to range
damm_tran <- function(
  x,hyper=list(
    KmO2_Rh=c(0.01,0.21,0),
    alphaVmaxSx=c(1666666.66666667,166666666.666667,1),
    EaVmaxSx=c(30,100,0),
    KMSx=c(8.00E-08,8.00E-06,1),
 #CH4  
    alphaVmaxCH4prod=c(1722222222222.22,172222222222222,1),
    EaVmaxCH4prod=c(15,75,0),
    KMSx_CH4prod=c(1666666.66666667,166666666.666667,1),
    alphaVmaxCH4ox=c(8.32660531276481,832.660531276481,1),
    EaVmaxCH4ox=c(10,50,0),
    KMCH4ox=c(5,50,0),
    KmO2_CH4ox=c(0.01,0.21,1),
    kl_CH4prod=c(0.001,0.01,1),
  #for N2O  
    alphaVmaxN2Oprod_nitrif=c(10,10000,1),
    EaVmaxN2Oprod_nitrif=c(25,200,1),
    KM_N2Oprod_nitrif=c(2.22,222.22,1),
    KMO2_N2Oprod_nitrif=c(14.33,1433.3,1),
    alphaVmaxN2Oprod_denitrif=c(10,10000,1),
    EaVmaxN2Oprod_denitrif=c(2.6,260,1),
    KM_N2Oprod_denitrif=c(2.22,222.22,1),
    alphaVmaxN2Ored=c(.0015,150000,1),
    EaVmaxN2Ored=c(10,100,0),
    KMN2Ored=c(1.00E-06,16.00,1),
    kl_N2Oprod=c(1.43,143.0,1),
    kl_N2Ored=c(0.72,72,1),
    KMSx_denitrif=c(10,100000,1)
    )
)
{
  ## x: unrestricted parameter on real line
  ## hyper: the range and whether or not to apply log transformation 
  ## value: parameter on restricted range
  names(x) <- NULL
  lKmO2_Rh <- x[1]
  lalphaVmaxSx <- x[2]
  lEaVmaxSx <- x[3]
  lKMSx <- x[4]
  lalphaVmaxCH4prod <- x[5]
  lEaVmaxCH4prod <- x[6]
  lKMSx_CH4prod <- x[7]
  lalphaVmaxCH4ox <- x[8]
  lEaVmaxCH4ox <- x[9]
  lKMCH4ox <- x[10]
  lKmO2_CH4ox <- x[11]
  lkl_CH4prod <- x[12]
  
  lalphaVmaxN2Oprod_nitrif <- x[13]
  lEaVmaxN2Oprod_nitrif <- x[14]
  lKM_N2Oprod_nitrif <- x[15]
  lKMO2_N2Oprod_nitrif <- x[16]
  lalphaVmaxN2Oprod_denitrif <- x[17]
  lEaVmaxN2Oprod_denitrif <- x[18]
  lKM_N2Oprod_denitrif <- x[19]
  lalphaVmaxN2Ored <- x[20]
  lEaVmaxN2Ored <- x[21]
  lKMN2Ored <- x[22]
  lkl_N2Oprod <- x[23]
  lkl_N2Ored <- x[24]
  lKMSx_denitrif <- x[25]
  ## KmO2_Rh
  tmp <- hyper[["KmO2_Rh"]]
  KmO2_Rh <- expit_tran(lKmO2_Rh,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxSx
  tmp <- hyper[["alphaVmaxSx"]]
  alphaVmaxSx <- expit_tran(lalphaVmaxSx,tmp[1],tmp[2],log=as.logical(tmp[3])) 
  ## EaVmaxSx
  tmp <- hyper[["EaVmaxSx"]]
  EaVmaxSx <- expit_tran(lEaVmaxSx,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KMSx
  tmp <- hyper[["KMSx"]]
  KMSx <- expit_tran(lKMSx,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxCH4prod
  tmp <- hyper[["alphaVmaxCH4prod"]]
  alphaVmaxCH4prod <- expit_tran(lalphaVmaxCH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## EaVmaxCH4prod
  tmp <- hyper[["EaVmaxCH4prod"]]
  EaVmaxCH4prod <- expit_tran(lEaVmaxCH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KMSx_CH4prod
  tmp <- hyper[["KMSx_CH4prod"]]
  KMSx_CH4prod <- expit_tran(lKMSx_CH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxCH4ox
  tmp <- hyper[["alphaVmaxCH4ox"]]
  alphaVmaxCH4ox <- expit_tran(lalphaVmaxCH4ox,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## EaVmaxCH4ox
  tmp <- hyper[["EaVmaxCH4ox"]]
  EaVmaxCH4ox <- expit_tran(lEaVmaxCH4ox,tmp[1],tmp[2],log=as.logical(tmp[3])) 
  ## KMCH4ox
  tmp <- hyper[["KMCH4ox"]]
  KMCH4ox <- expit_tran(lKMCH4ox,tmp[1],tmp[2],log=as.logical(tmp[3])) 
  ## KmO2_CH4ox
  tmp <- hyper[["KmO2_CH4ox"]]
  KmO2_CH4ox <- expit_tran(lKmO2_CH4ox,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## kl_CH4prod
  tmp <- hyper[["kl_CH4prod"]]
  kl_CH4prod <- expit_tran(lkl_CH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  
  ## alphaVmaxN2Oprod_nitrif
  tmp <- hyper[["alphaVmaxN2Oprod_nitrif"]]
  alphaVmaxN2Oprod_nitrif <- expit_tran(lalphaVmaxN2Oprod_nitrif,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## EaVmaxN2Oprod_nitrif
  tmp <- hyper[["EaVmaxN2Oprod_nitrif"]]
  EaVmaxN2Oprod_nitrif <- expit_tran(lEaVmaxN2Oprod_nitrif,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KM_N2Oprod_nitrif
  tmp <- hyper[["KM_N2Oprod_nitrif"]]
  KM_N2Oprod_nitrif <- expit_tran(lKM_N2Oprod_nitrif,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KMO2_N2Oprod_nitrif
  tmp <- hyper[["KMO2_N2Oprod_nitrif"]]
  KMO2_N2Oprod_nitrif <- expit_tran(lKMO2_N2Oprod_nitrif,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxN2Oprod_denitrif
  tmp <- hyper[["alphaVmaxN2Oprod_denitrif"]]
  alphaVmaxN2Oprod_denitrif <- expit_tran(lalphaVmaxN2Oprod_denitrif,tmp[1],tmp[2],log=as.logical(tmp[3])) 
  ## EaVmaxN2Oprod_denitrif
  tmp <- hyper[["EaVmaxN2Oprod_denitrif"]]
  EaVmaxN2Oprod_denitrif <- expit_tran(lEaVmaxN2Oprod_denitrif,tmp[1],tmp[2],log=as.logical(tmp[3])) 
  ## KM_N2Oprod_denitrif
  tmp <- hyper[["KM_N2Oprod_denitrif"]]
  KM_N2Oprod_denitrif <- expit_tran(lKM_N2Oprod_denitrif,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxN2Ored
  tmp <- hyper[["alphaVmaxN2Ored"]]
  alphaVmaxN2Ored <- expit_tran(lalphaVmaxN2Ored,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## EaVmaxN2Ored
  tmp <- hyper[["EaVmaxN2Ored"]]
  EaVmaxN2Ored <- expit_tran(lEaVmaxN2Ored,tmp[1],tmp[2],log=as.logical(tmp[3]))  
  ## KMN2Ored
  tmp <- hyper[["KMN2Ored"]]
  KMN2Ored <- expit_tran(lKMN2Ored,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## kl_N2Oprod
  tmp <- hyper[["kl_N2Oprod"]]
  kl_N2Oprod <- expit_tran(lkl_N2Oprod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## kl_N2Ored
  tmp <- hyper[["kl_N2Ored"]]
  kl_N2Ored <- expit_tran(lkl_N2Ored,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KMSx_denitrif
  tmp <- hyper[["KMSx_denitrif"]]
  KMSx_denitrif <- expit_tran(lKMSx_denitrif,tmp[1],tmp[2],log=as.logical(tmp[3]))
  
  r <- c(KmO2_Rh=KmO2_Rh,alphaVmaxSx=alphaVmaxSx,EaVmaxSx=EaVmaxSx,KMSx=KMSx,
         alphaVmaxCH4prod=alphaVmaxCH4prod,EaVmaxCH4prod=EaVmaxCH4prod,
         KMSx_CH4prod=KMSx_CH4prod,alphaVmaxCH4ox=alphaVmaxCH4ox,
         EaVmaxCH4ox=EaVmaxCH4ox,KMCH4ox=KMCH4ox,
         KmO2_CH4ox=KmO2_CH4ox,kl_CH4prod=kl_CH4prod,
         alphaVmaxN2Oprod_nitrif=alphaVmaxN2Oprod_nitrif, EaVmaxN2Oprod_nitrif=EaVmaxN2Oprod_nitrif,  
         KM_N2Oprod_nitrif=KM_N2Oprod_nitrif,KMO2_N2Oprod_nitrif=KMO2_N2Oprod_nitrif,
         alphaVmaxN2Oprod_denitrif=alphaVmaxN2Oprod_denitrif, EaVmaxN2Oprod_denitrif=EaVmaxN2Oprod_denitrif,
         KM_N2Oprod_denitrif=KM_N2Oprod_denitrif,alphaVmaxN2Ored=alphaVmaxN2Ored, 
         EaVmaxN2Ored=EaVmaxN2Ored, KMN2Ored=KMN2Ored, 
         kl_N2Oprod=kl_N2Oprod,kl_N2Ored=kl_N2Ored,KMSx_denitrif=KMSx_denitrif)
  #names(r) <- names(x)
  r
}

## helper function to back transform restricted parameters to real
damm_btran <- function(
  z,hyper=list(
    KmO2_Rh=c(0.01,0.21,0),
    alphaVmaxSx=c(1666666.66666667,166666666.666667,1),
    EaVmaxSx=c(30,100,0),
    KMSx=c(8.00E-08,8.00E-06,1),
    alphaVmaxCH4prod=c(1722222222222.22,172222222222222,1),
    EaVmaxCH4prod=c(15,75,0),
    KMSx_CH4prod=c(1666666.66666667,166666666.666667,1),
    alphaVmaxCH4ox=c(8.32660531276481,832.660531276481,1),
    EaVmaxCH4ox=c(10,50,0),KMCH4ox=c(5,50,0),
    KmO2_CH4ox=c(0.01,0.21,1),
    kl_CH4prod=c(0.001,0.01,1),
    alphaVmaxN2Oprod_nitrif=c(10,10000,1),
    EaVmaxN2Oprod_nitrif=c(25,200,1),
    KM_N2Oprod_nitrif=c(2.22,222.22,1),
    KMO2_N2Oprod_nitrif=c(14.33,1433.3,1),
    alphaVmaxN2Oprod_denitrif=c(10,10000,1),
    EaVmaxN2Oprod_denitrif=c(2.6,260,1),
    KM_N2Oprod_denitrif=c(2.22,222.22,1),
    alphaVmaxN2Ored=c(.0015,150000,1),
    EaVmaxN2Ored=c(10,100,0),
    KMN2Ored=c(1.00E-06,16.00,1),
    kl_N2Oprod=c(1.43,143.0,1),
    kl_N2Ored=c(0.72,72,1),
    KMSx_denitrif=c(10,100000,1)
  )
)
{
  ## z: parameter on uniform range
  ## value: parameter on real line
  names(z) <- NULL
  
  KmO2_Rh <- z[1]
  alphaVmaxSx <- z[2]
  EaVmaxSx <- z[3]
  KMSx <- z[4]
  alphaVmaxCH4prod <- z[5]
  EaVmaxCH4prod <- z[6]
  KMSx_CH4prod <- z[7]
  alphaVmaxCH4ox <- z[8]
  EaVmaxCH4ox <- z[9]
  KMCH4ox <- z[10]
  KmO2_CH4ox <- z[11]
  kl_CH4prod <- z[12]
  
  alphaVmaxN2Oprod_nitrif <- z[13]
  EaVmaxN2Oprod_nitrif <- z[14]
  KM_N2Oprod_nitrif <- z[15]
  KMO2_N2Oprod_nitrif <- z[16]
  alphaVmaxN2Oprod_denitrif <- z[17]
  EaVmaxN2Oprod_denitrif <- z[18]
  KM_N2Oprod_denitrif <- z[19]
  alphaVmaxN2Ored <- z[20]
  EaVmaxN2Ored <- z[21]
  KMN2Ored <- z[22]
  kl_N2Oprod <- z[23]
  kl_N2Ored <- z[24]
  KMSx_denitrif <- z[25]
  
  ## KmO2_Rh
  tmp <- hyper[["KmO2_Rh"]]
  lKmO2_Rh <- logit_tran(KmO2_Rh,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxSx
  tmp <- hyper[["alphaVmaxSx"]]
  lalphaVmaxSx=logit_tran(alphaVmaxSx,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## EaVmaxSx
  tmp <- hyper[["EaVmaxSx"]]
  lEaVmaxSx=logit_tran(EaVmaxSx,tmp[1],tmp[2],log=as.logical(tmp[3])) 
  ## KMSx
  tmp <- hyper[["KMSx"]]
  lKMSx=logit_tran(KMSx,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxCH4prod
  tmp <- hyper[["alphaVmaxCH4prod"]]
  lalphaVmaxCH4prod <- logit_tran(alphaVmaxCH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## EaVmaxCH4prod
  tmp <- hyper[["EaVmaxCH4prod"]]
  lEaVmaxCH4prod=logit_tran(EaVmaxCH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KMSx_CH4prod
  tmp <- hyper[["KMSx_CH4prod"]]
  lKMSx_CH4prod <- logit_tran(KMSx_CH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxCH4ox
  tmp <- hyper[["alphaVmaxCH4ox"]]
  lalphaVmaxCH4ox <- logit_tran(alphaVmaxCH4ox,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## EaVmaxCH4ox
  tmp <- hyper[["EaVmaxCH4ox"]]
  lEaVmaxCH4ox=logit_tran(EaVmaxCH4ox,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KMCH4ox
  tmp <- hyper[["KMCH4ox"]]
  lKMCH4ox=logit_tran(KMCH4ox,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KmO2_CH4ox
  tmp <- hyper[["KmO2_CH4ox"]]
  lKmO2_CH4ox <- logit_tran(KmO2_CH4ox,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## kl_CH4prod
  tmp <- hyper[["kl_CH4prod"]]
  lkl_CH4prod <- logit_tran(kl_CH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  
  ## alphaVmaxN2Oprod_nitrif
  tmp <- hyper[["alphaVmaxN2Oprod_nitrif"]]
  lalphaVmaxN2Oprod_nitrif <- logit_tran(alphaVmaxN2Oprod_nitrif,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## EaVmaxN2Oprod_nitrif
  tmp <- hyper[["EaVmaxN2Oprod_nitrif"]]
  lEaVmaxN2Oprod_nitrif <- logit_tran(EaVmaxN2Oprod_nitrif,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KM_N2Oprod_nitrif
  tmp <- hyper[["KM_N2Oprod_nitrif"]]
  lKM_N2Oprod_nitrif <- logit_tran(KM_N2Oprod_nitrif,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KMO2_N2Oprod_nitrif
  tmp <- hyper[["KMO2_N2Oprod_nitrif"]]
  lKMO2_N2Oprod_nitrif <- logit_tran(KMO2_N2Oprod_nitrif,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxN2Oprod_denitrif
  tmp <- hyper[["alphaVmaxN2Oprod_denitrif"]]
  lalphaVmaxN2Oprod_denitrif <- logit_tran(alphaVmaxN2Oprod_denitrif,tmp[1],tmp[2],log=as.logical(tmp[3])) 
  ## EaVmaxN2Oprod_denitrif
  tmp <- hyper[["EaVmaxN2Oprod_denitrif"]]
  lEaVmaxN2Oprod_denitrif <- logit_tran(EaVmaxN2Oprod_denitrif,tmp[1],tmp[2],log=as.logical(tmp[3])) 
  ## KM_N2Oprod_denitrif
  tmp <- hyper[["KM_N2Oprod_denitrif"]]
  lKM_N2Oprod_denitrif <- logit_tran(KM_N2Oprod_denitrif,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxN2Ored
  tmp <- hyper[["alphaVmaxN2Ored"]]
  lalphaVmaxN2Ored <- logit_tran(alphaVmaxN2Ored,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## EaVmaxN2Ored
  tmp <- hyper[["EaVmaxN2Ored"]]
  lEaVmaxN2Ored <- logit_tran(EaVmaxN2Ored,tmp[1],tmp[2],log=as.logical(tmp[3]))  
  ## KMN2Ored
  tmp <- hyper[["KMN2Ored"]]
  lKMN2Ored <- logit_tran(KMN2Ored,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## kl_N2Oprod
  tmp <- hyper[["kl_N2Oprod"]]
  lkl_N2Oprod <- logit_tran(kl_N2Oprod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## kl_N2Ored
  tmp <- hyper[["kl_N2Ored"]]
  lkl_N2Ored <- logit_tran(kl_N2Ored,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KMSx_denitrif
  tmp <- hyper[["KMSx_denitrif"]]
  lKMSx_denitrif <- logit_tran(KMSx_denitrif,tmp[1],tmp[2],log=as.logical(tmp[3]))
  
  r <- c(lKmO2_Rh=lKmO2_Rh,
         lalphaVmaxSx=lalphaVmaxSx,lEaVmaxSx=lEaVmaxSx,lKMSx=lKMSx,
         lalphaVmaxCH4prod=lalphaVmaxCH4prod,lEaVmaxCH4prod=lEaVmaxCH4prod,
         lKMSx_CH4prod=lKMSx_CH4prod,lalphaVmaxCH4ox=lalphaVmaxCH4ox,
         lEaVmaxCH4ox=lEaVmaxCH4ox,lKMCH4ox=lKMCH4ox,
         lKmO2_CH4ox=lKmO2_CH4ox,lkl_CH4prod=lkl_CH4prod,
         lalphaVmaxN2Oprod_nitrif=lalphaVmaxN2Oprod_nitrif,lEaVmaxN2Oprod_nitrif=lEaVmaxN2Oprod_nitrif,  
         lKM_N2Oprod_nitrif=lKM_N2Oprod_nitrif,lKMO2_N2Oprod_nitrif=lKMO2_N2Oprod_nitrif,
         lalphaVmaxN2Oprod_denitrif=lalphaVmaxN2Oprod_denitrif,
         lEaVmaxN2Oprod_denitrif=lEaVmaxN2Oprod_denitrif,
         lKM_N2Oprod_denitrif=lKM_N2Oprod_denitrif,lalphaVmaxN2Ored=lalphaVmaxN2Ored, 
         lEaVmaxN2Ored=lEaVmaxN2Ored, lKMN2Ored=lKMN2Ored, 
         lkl_N2Oprod=lkl_N2Oprod,lkl_N2Ored=lkl_N2Ored,lKMSx_denitrif=lKMSx_denitrif)
  #names(r) <- names(z)
  r
}

tran_test_08082017 <- function()
{
  #setwd("C:/Users/dsihi/Google Drive/sihidong")
  rm(list=ls())
  source("damm_utils_3.R")
  hyper=list(
    KmO2_Rh=c(0.01,0.21,0),
    alphaVmaxSx=c(1666666.66666667,166666666.666667,1),
    EaVmaxSx=c(30,100,0),
    KMSx=c(8.00E-08,8.00E-06,1),
    alphaVmaxCH4prod=c(1722222222222.22,172222222222222,1),
    EaVmaxCH4prod=c(15,75,0),
    KMSx_CH4prod=c(1666666.66666667,166666666.666667,1),
    alphaVmaxCH4ox=c(8.32660531276481,832.660531276481,1),
    EaVmaxCH4ox=c(10,50,0),KMCH4ox=c(5,50,0),
    KmO2_CH4ox=c(0.01,0.21,1),
    kl_CH4prod=c(0.001,0.01,1),
    
    alphaVmaxN2Oprod_nitrif=c(10,10000,1),
    EaVmaxN2Oprod_nitrif=c(25,200,1),
    KM_N2Oprod_nitrif=c(2.22,222.22,1),
    KMO2_N2Oprod_nitrif=c(14.33,1433.3,1),
    alphaVmaxN2Oprod_denitrif=c(10,10000,1),
    EaVmaxN2Oprod_denitrif=c(2.6,260,1),
    KM_N2Oprod_denitrif=c(2.22,222.22,1),
    alphaVmaxN2Ored=c(.0015,150000,1),
    EaVmaxN2Ored=c(10,100,0),
    KMN2Ored=c(1.00E-06,16.00,1),
    kl_N2Oprod=c(1.43,143.0,1),
    kl_N2Ored=c(0.72,72,1),
    KMSx_denitrif=c(10,100000,1)
  )
  x1=c(
    KmO2_Rh=1.417843,
    alphaVmaxSx=15.1561754430518,
    EaVmaxSx=0.224663729,
    KMSx=-15.46972856,
    alphaVmaxCH4prod=15.56470445,
    EaVmaxCH4prod=1.086290931,
    KMSx_CH4prod=-15.46972856,
    alphaVmaxCH4ox=5.130463522,
    EaVmaxCH4ox=-0.025674405,
    KMCH4ox=-2.77212177,
    KmO2_CH4ox=0,
    kl_CH4prod=0,
    
    alphaVmaxN2Oprod_nitrif=100, 
    EaVmaxN2Oprod_nitrif=50,  
    KM_N2Oprod_nitrif=22.22,
    KMO2_N2Oprod_nitrif=143.3,
    alphaVmaxN2Oprod_denitrif=100, 
    EaVmaxN2Oprod_denitrif=50,  
    KM_N2Oprod_denitrif=26,
    alphaVmaxN2Ored=1500, 
    EaVmaxN2Ored=50,   
    KMN2Ored=0.16, 
    kl_N2Oprod=14.3,
    kl_N2Ored=7.2,
    KMSx_denitrif=1000
  )
  
  damm_tran(x1)
  
  # KmO2_Rh      alphaVmaxSx         EaVmaxSx             KMSx alphaVmaxCH4prod 
  # 1.710000e-01     1.666665e+08     6.891516e+01     8.000007e-08     1.722221e+14 
  # EaVmaxCH4prod     KMSx_CH4prod   alphaVmaxCH4ox      EaVmaxCH4ox          KMCH4ox 
  # 5.986096e+01     1.666668e+06     8.104195e+02     2.974327e+01     7.648222e+00 
  # KmO2_CH4ox       kl_CH4prod 
  # 4.582576e-02     3.162278e-03 
  x0 <-   c(
    lKmO2_Rh=logit_tran(0.1,0.01,0.21),
    lalphaVmaxSx=logit_tran(16666666.6666667,hyper[["alphaVmaxSx"]][1],hyper[["alphaVmaxSx"]][2],log=TRUE), 
    lEaVmaxSx=logit_tran(68,hyper[["EaVmaxSx"]][1],hyper[["EaVmaxSx"]][2]),
    lKMSx=logit_tran(0.0000008,hyper[["KMSx"]][1],hyper[["KMSx"]][2],log=TRUE),
    lalphaVmaxCH4prod=logit_tran((6.2*(10^16)*(1/3600)),hyper[["alphaVmaxCH4prod"]][1],hyper[["alphaVmaxCH4prod"]][2],log=TRUE), 
    lEaVmaxCH4prod=logit_tran(50,hyper[["EaVmaxCH4prod"]][1],hyper[["EaVmaxCH4prod"]][2]),
    lKMSx_CH4prod=logit_tran(2*16666666.6666667,1666666.66666667,166666666.666667,log=TRUE),
    lalphaVmaxCH4ox=logit_tran(299757.791259533*(1/3600),hyper[["alphaVmaxCH4ox"]][1],hyper[["alphaVmaxCH4ox"]][2],log=TRUE), 
    lEaVmaxCH4ox=logit_tran(30,hyper[["EaVmaxCH4ox"]][1],hyper[["EaVmaxCH4ox"]][2]),
    lKMCH4ox=logit_tran(7,hyper[["KMCH4ox"]][1],hyper[["KMCH4ox"]][2]),
    lKmO2_CH4ox=logit_tran(0.10,0.01,0.21),
    lkl_CH4prod=logit_tran(0.005,0.001,0.01,log=TRUE)
  )
  options(digits=2,warn=1)
  (z <- damm_tran(x0))
  # KmO2_Rh      alphaVmaxSx         EaVmaxSx             KMSx alphaVmaxCH4prod 
  # 1.0e-01          1.7e+07          6.8e+01          8.0e-07          1.7e+13 
  # EaVmaxCH4prod     KMSx_CH4prod   alphaVmaxCH4ox      EaVmaxCH4ox          KMCH4ox 
  # 5.0e+01          3.3e+07          8.3e+01          3.0e+01          7.0e+00 
  # KmO2_CH4ox       kl_CH4prod 
  # 3.9e-02          5.0e-03 
  (x1 <- damm_btran(z))
  x0-x1
  # lKmO2_Rh      lalphaVmaxSx         lEaVmaxSx             lKMSx 
  # 0.0e+00           0.0e+00           0.0e+00           1.6e-15 
  # lalphaVmaxCH4prod    lEaVmaxCH4prod     lKMSx_CH4prod   lalphaVmaxCH4ox 
  # 0.0e+00           0.0e+00           0.0e+00           0.0e+00 
  # lEaVmaxCH4ox          lKMCH4ox       lKmO2_CH4ox       lkl_CH4prod 
  # 0.0e+00           0.0e+00           1.4e-16           0.0e+00 
  
  z <- c(KmO2_Rh=0.10,
         alphaVmaxSx=2666666.66666667,
         EaVmaxSx=70,
         KMSx=5.00E-06,
         alphaVmaxCH4prod=(6.2*(10^16)*(1/3600)),
         EaVmaxCH4prod=25,
         KMSx_CH4prod=2*16666666.6666667,
         alphaVmaxCH4ox=83.2660531276481,
         EaVmaxCH4ox=30,
         KMCH4ox=7,
         KmO2_CH4ox=0.10,
         kl_CH4prod=0.005)
  (x <- damm_btran(z))
  (z2 <- damm_tran(x))
  z2-z
  # KmO2_Rh      alphaVmaxSx         EaVmaxSx             KMSx alphaVmaxCH4prod 
  # 0.0e+00         -2.3e-09          0.0e+00         -1.7e-21         -2.5e-02 
  # EaVmaxCH4prod     KMSx_CH4prod   alphaVmaxCH4ox      EaVmaxCH4ox          KMCH4ox 
  # 0.0e+00         -4.8e-08          0.0e+00          0.0e+00          0.0e+00 
  # KmO2_CH4ox       kl_CH4prod 
  # 1.4e-17          1.7e-18
  
}