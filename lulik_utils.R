## lulik_utils.R
  
## utility function to run non-linear least square 
## and the conventional Nelder-Mead algorithm using optim

## most arugments are defined in damm_xx.R and lulik_xx.R
## the first two arugments have been modified to allow piece-wise optimization
## par: the interested paramemter
## z  : the ancillary parameter
## for rh, CH4 N2O  updating was not performed
## value: the residual as defined by minpack.lm
#lulik.rh return residuals as vector
lulik.rh <- function(
  par,z,formula,data,nlevels,
  R, O2airfrac, CH4airfrac, N2Oairfrac,
  BD, PD,  Porosity,
  p, RQ,  Depth_O2, Depth_CH4,Depth_N2O,
  microsite,total.microsite,time.step,
  flagO2=TRUE,  flagCH4=FALSE, flagN2O=FALSE, ## whether run O2 and CH4 N2O to steady states
  iO2.lst=NULL,iCH4.lst=NULL,iN2O.lst=NULL,
  reltol=sqrt(.Machine$double.eps),verbose=0,
  sxtot.range=c(0.01,0.15), nquantile=10, mult=1/4,
  soilm.cv=20, soilm.range=c(70,130), soilm.nlevel=10, lognorm=TRUE,
  nu=2,R0=diag(2),
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
  ))
{
  r <- lulik(
    x=c(par,z),formula=formula,data=data,nlevels=nlevels,
    R=R, O2airfrac=O2airfrac, CH4airfrac = CH4airfrac,N2Oairfrac = N2Oairfrac,
    BD=BD, PD=PD,  Porosity=Porosity,
    p=p, RQ=RQ,  Depth_O2 = Depth_O2, Depth_CH4 = Depth_CH4,Depth_N2O = Depth_N2O,
    microsite=microsite,total.microsite=total.microsite,time.step=time.step,
    flagO2=flagO2,  flagCH4=flagCH4, flagN2O=flagN2O, 
    iO2.lst=iO2.lst,iCH4.lst=iCH4.lst,iN2O.lst=iN2O.lst,
    reltol=reltol,verbose=verbose,
    sxtot.range=sxtot.range, nquantile=nquantile, mult=mult,
    soilm.cv=soilm.cv, soilm.range=soilm.range, soilm.nlevel=soilm.nlevel, 
    lognorm=lognorm,
    nu=nu,R0=R0,
    hyper=hyper)
  
  resid <- r$data[,3] - r$fit_Rh
  resid[!is.na(resid)]
}
#.optim return square error(sum of residual^2)
lulik.rh.optim <- function(
  par,z,formula,data,nlevels,
  R, O2airfrac, CH4airfrac, N2Oairfrac,BD, PD,  Porosity,
  p, RQ,  Depth_O2, Depth_CH4,Depth_N2O,
  microsite,total.microsite,time.step,
  flagO2=TRUE,  flagCH4=FALSE, flagN2O=FALSE, ## whether run O2 and CH4 N2O to steady states
  iO2.lst=NULL,iCH4.lst=NULL,iN2O.lst=NULL,
  reltol=sqrt(.Machine$double.eps),verbose=0,
  sxtot.range=c(0.01,0.15), nquantile=10, mult=1/4,
  soilm.cv=20, soilm.range=c(70,130), soilm.nlevel=10, lognorm=TRUE,
  nu=2,R0=diag(2),
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
  ))
{
  r <- lulik(
    x=c(par,z),formula=formula,data=data,nlevels=nlevels,
    R=R, O2airfrac=O2airfrac, CH4airfrac = CH4airfrac,N2Oairfrac = N2Oairfrac,
    BD=BD, PD=PD,  Porosity=Porosity,
    p=p, RQ=RQ,  Depth_O2=Depth_O2, Depth_CH4 = Depth_CH4,Depth_N2O = Depth_N2O, 
    microsite=microsite,total.microsite=total.microsite,time.step=time.step,
    flagO2=flagO2,  flagCH4=flagCH4,  flagN2O=flagN2O, 
    iO2.lst=iO2.lst,iCH4.lst=iCH4.lst,iN2O.lst=iN2O.lst,
    reltol=reltol,verbose=verbose,
    sxtot.range=sxtot.range, nquantile=nquantile, mult=mult,
    soilm.cv=soilm.cv, soilm.range=soilm.range, soilm.nlevel=soilm.nlevel, 
    lognorm=lognorm,
    nu=nu,R0=R0,
    hyper=hyper)
  
  resid <- r$data[,3] - r$fit_Rh
  resid <- resid[!is.na(resid)]
  as.vector(crossprod(resid))
}

lulik.ch4.n2o <- function(
  par,z,formula,data,nlevels,
  R, O2airfrac, CH4airfrac, N2Oairfrac,BD, PD, Porosity,
  p, RQ,  Depth_O2, Depth_CH4,Depth_N2O,
  microsite,total.microsite,time.step,
  flagO2=TRUE,  flagCH4=TRUE, flagN2O=TRUE,## whether run O2 and CH4 to steady states
  iO2.lst=NULL,iCH4.lst=NULL,iN2O.lst=NULL,
  reltol=sqrt(.Machine$double.eps),verbose=0,
  sxtot.range=c(0.01,0.15), nquantile=10, mult=1/4,
  soilm.cv=20, soilm.range=c(70,130), soilm.nlevel=10, lognorm=TRUE,
  nu=2,R0=diag(2),
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
  ))
{
  r <- lulik(
    x=c(z,par),formula=formula,data=data,nlevels=nlevels,
    R=R, O2airfrac=O2airfrac, CH4airfrac = CH4airfrac, N2Oairfrac = N2Oairfrac, 
    BD=BD, PD=PD,  Porosity=Porosity,
    p=p, RQ=RQ,  Depth_O2=Depth_O2, Depth_CH4 = Depth_CH4,Depth_N2O = Depth_N2O,
    microsite=microsite,total.microsite=total.microsite,time.step=time.step,
    flagO2=flagO2,  flagCH4=flagCH4,  flagN2O=flagN2O, 
    iO2.lst=iO2.lst,iCH4.lst=iCH4.lst,iN2O.lst=iN2O.lst,
    reltol=reltol,verbose=verbose,
    sxtot.range=sxtot.range, nquantile=nquantile, mult=mult,
    soilm.cv=soilm.cv, soilm.range=soilm.range, soilm.nlevel=soilm.nlevel, 
    lognorm=lognorm,
    nu=nu,R0=R0,
    hyper=hyper)
  
  resid <- c((r$data[,4] - r$fit_CH4),(r$data[,5] - r$fit_N2O))
  resid[!is.na(resid)]
  
    # resid <- r$data[,5] - r$fit_N2O
    # resid[!is.na(resid)]
}

lulik.ch4.n2o.optim <- function(
  par,z,formula,data,nlevels,
  R, O2airfrac, CH4airfrac, N2Oairfrac,BD, PD,  Porosity,
  p, RQ,  Depth_O2, Depth_CH4,Depth_N2O,
  microsite,total.microsite,time.step,
  flagO2=TRUE,  flagCH4=TRUE, flagN2O=TRUE,## whether run O2 and CH4 to steady states
  iO2.lst=NULL,iCH4.lst=NULL,iN2O.lst=NULL,
  reltol=sqrt(.Machine$double.eps),verbose=0,
  sxtot.range=c(0.01,0.15), nquantile=10, mult=1/4,
  soilm.cv=20, soilm.range=c(70,130), soilm.nlevel=10, lognorm=TRUE,
  nu=2,R0=diag(2),
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
  ))
{
  r <- lulik(
    x=c(z,par),formula=formula,data=data,nlevels=nlevels,
    R=R, O2airfrac=O2airfrac, CH4airfrac = CH4airfrac,N2Oairfrac = N2Oairfrac,
    BD=BD, PD=PD,  Porosity=Porosity,
    p=p, RQ=RQ,  Depth_O2=Depth_O2, Depth_CH4 = Depth_CH4, Depth_N2O = Depth_N2O,
    microsite=microsite,total.microsite=total.microsite,time.step=time.step,
    flagO2=flagO2,  flagCH4=flagCH4,  flagN2O=flagN2O, 
    iO2.lst=iO2.lst,iCH4.lst=iCH4.lst,iN2O.lst=iN2O.lst,
    reltol=reltol,verbose=verbose,
    sxtot.range=sxtot.range, nquantile=nquantile, mult=mult,
    soilm.cv=soilm.cv, soilm.range=soilm.range, soilm.nlevel=soilm.nlevel, 
    lognorm=lognorm,
    nu=nu,R0=R0,
    hyper=hyper)
  
  resid <- c((r$data[,4] - r$fit_CH4),(r$data[,5] - r$fit_N2O))
  resid <- resid[!is.na(resid)]
  as.vector(crossprod(resid))
}


#not run, need to add N2O parameter and N2O
test <- function()
{
  rm(list=ls())
  ## Session > Work Dir > Source File Loc
  trench <- read.csv("2015data.csv",head=T)
  # 'data.frame':	111 obs. of  5 variables:
  #   $ Data..point: int  1 2 3 4 5 6 7 8 9 10 ...
  # $ SoilM      : num  28 27.8 32.3 30.1 29 ...
  # $ SoilT      : num  12.8 13.3 14.3 13.6 12.8 ...
  # $ Rh         : num  145 145 158 131 113 ...
  # $ CH4        : num  -0.0407 -0.0428 -0.0356 -0.0375 -0.0349 .
  ## Initial values to start estimation
  z0 <-   c(
    KmO2_Rh=245,alphaVmaxSx=230000000000, EaVmaxSx=72,
    KMSx=1000,
    alphaVmaxCH4prod=3000000,  EaVmaxCH4prod=100,
    KMSx_CH4prod=350,
    alphaVmaxCH4ox=7, EaVmaxCH4ox=30,
    KMCH4ox=1.00E-02,  KmO2_CH4ox=43,kl_CH4prod=3
  )
  
  
  hyper=list(
    KmO2_Rh=c(2.45,24500,1),
    alphaVmaxSx=c(2300000000,23000000000000,1),
    EaVmaxSx=c(7.2,7200,1),
    KMSx=c(10,100000,1),
    alphaVmaxCH4prod=c(30000,300000000,1),
    EaVmaxCH4prod=c(0.3,3000,1),
    KMSx_CH4prod=c(3.5,35000,1),
    alphaVmaxCH4ox=c(0.07,700,1),
    EaVmaxCH4ox=c(10,50,0),KMCH4ox=c(1.00E-04,1.00,1),
    KmO2_CH4ox=c(0.043,4300,1),
    kl_CH4prod=c(0.03,300,1)
  )
  
  source("damm_11a.R")
  source("damm_utils_3.R")
  x0 <-   damm_btran(z0,hyper)

  source("lulik_4b.R")
  source("lulik_utils_2b.R")
  r0 <- lulik.rh(
    par=x0[1:4],z=x0[5:12],formula=~SoilM+SoilT+Rh+CH4+N2O+NH4+NO3,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8,
    BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10,
    microsite = 0,total.microsite = 10000, time.step = 300
  )
  hist(r0)
  
  r1 <- lulik.ch4(
    par=x0[5:12],z=x0[1:4],formula=~SoilM+SoilT+Rh+CH4+N2O+NH4+NO3,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8,
    BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10,
    microsite = 0,total.microsite = 10000, time.step = 300
  )
  hist(r1)
}