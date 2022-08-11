## lulik.R

## log likelihood function given data, 
## marginalized over the precision matrix
## ver 4) allow missing inputs, no longer need marginal likelihood calculation
#source("mlik_1.R")
require(doParallel)
## arguments:
## x : parameter on the un-restricted scale - real line
## formula: a one sided formula defining the observed soil moisture
##         soil temperature, CO2 and methane in that order
## data : input data frame
## hyper: defining the range and whether log-transformation is needed
## nlevels: number of levels for bivariate distribution (check consistency)
## iO2.lst: a list of initial O2 values (default 0.209?)
## iCH4.lst: a list of initial CH4 values (default 1.8?)
## iN2O.lst: a list of initial N2O values (default 1.8?)
## nu   : prior degree of freedom for 2x2 precision (default=2)
## R0   : prior concentration matrix for 2x2 precision (default identity)
## value: a list with the following entries
##  value: log likelihood value to minimize
##  fit_Rh: model prediction of Rh
##  fit_CH4: model prediction of CH4
##  fit_N2O: model prediction of N2O
##  par: input parameter
##  tpar: input parameter on the uniform and restricted scale
##  flag: convergence flag
##  output: microsite and time specific steady state
##  O2.lst: output final O2 values
##  CH4.lst: output final CH4 values
##  N2O.lst: output final N2O values
lulik <- function(
  x,formula,data,nlevels,
  R, O2airfrac, CH4airfrac, N2Oairfrac,BD, PD,  Porosity,
  p, RQ,  Depth_O2, Depth_CH4,Depth_N2O,
  microsite,total.microsite,time.step,
  flagO2=TRUE,  flagCH4=TRUE, flagN2O=TRUE,## whether run O2 and CH4 N2Oto steady states
  iO2.lst=NULL,iCH4.lst=NULL,iN2O.lst=NULL,
  reltol=sqrt(.Machine$double.eps),verbose=0,
  sxtot.range=c(0.01,0.15), nquantile=10, mult=1/4,
  soilm.cv=20, soilm.range=c(70,200), soilm.nlevel=10, lognorm=TRUE,
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
    EaVmaxN2Ored=c(10,100,1),
    KMN2Ored=c(1.00E-06,16.00,1),
    kl_N2Oprod=c(1.43,143.0,1),
    kl_N2Ored=c(0.72,72,1),
    KMSx_denitrif=c(10,100000,1) 
  ))
{
  ts0 <- proc.time()
  ## Parameter Transformation
  # stopifnot(!any(is.na(x)))
  z <- damm_tran(x=x,hyper=hyper) ## back-transform to restricted range,z is parameter initial values
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
  
  ## extract observed data (alow NA in responses)
  #l.mat <- as.matrix(na.omit(data[,all.vars(formula)]))
  l.mat0 <- as.matrix(data[,all.vars(formula)])
  ii <- which(!is.na(l.mat0[,1])&!is.na(l.mat0[,2]))
  l.mat <- l.mat0[ii,]
  #stopifnot(all(!is.na(l.mat[,1])))
  #stopifnot(all(!is.na(l.mat[,2])))
  ## Note: missing values automatically dropped
  SoilM__ <- l.mat[,1]
  SoilT__ <- l.mat[,2]
  Rh__ <- l.mat[,3]
  CH4__ <- l.mat[,4]
  
  N2O__ <- l.mat[,5]
  
  NH4__ <- l.mat[,6] #THIS IS THE PLACE TO READ NH4 and NO3
  NO3__ <- l.mat[,7]
  
  
  #browser()
  ## steady state for each time point
  ## storage
  n <- nrow(l.mat)#each row is an observation
  meanRH <- rep(NA,n)
  meanCH4 <- rep(NA,n)
  meanN2O <- rep(NA,n)
  flag <- vector("list",n) ## convergence flags
  olst <- vector("list",n) ## microsite specific output
  O2.lst <- vector("list",n) ## final O2 values
  CH4.lst <- vector("list",n) ## final CH4 values
  N2O.lst <- vector("list",n) ## final N2O values
  big.lst <- foreach(i = 1:n) %dopar% {#this loop for each observation
  #for(i in 1:n){  
    if(verbose>1) cat(i,'\n')
    ## extract initial values
    if(!is.null(iO2.lst)){
      iO2 <- iO2.lst[[i]]
    }else{
      #iO2 <- rep(0.209,nlevels) # 11_a
      iO2 <- rep(O2airfrac,nlevels)
    }
    if(!is.null(iCH4.lst)){
      iCH4 <- iCH4.lst[[i]]
    }else{
      #iCH4 <- rep(1.8,nlevels) # 11_a
      iCH4 <- CH4airfrac
    }
    
    if(!is.null(iN2O.lst)){
      iN2O <- iN2O.lst[[i]]
    }else{
      #iN2O <- rep(0.328,nlevels) # 11_a
      iN2O <- N2Oairfrac
    }
    
    ## run steady state for time point t
    tmp <- damm(
      R=R, O2airfrac = O2airfrac, CH4airfrac=CH4airfrac,N2Oairfrac=N2Oairfrac,
      BD=BD, PD=PD,Porosity = Porosity,
      p=p, RQ=RQ,
      SoilM.o=SoilM__[i], SoilT = SoilT__[i],
      NH4= NH4__[i],NO3= NO3__[i],
      Depth_O2=Depth_O2, Depth_CH4=Depth_CH4,Depth_N2O=Depth_N2O,
      microsite = microsite, total.microsite=total.microsite,
      time.step = time.step, 
      iO2 = iO2, iCH4 = iCH4,iN2O = iN2O,
      flagO2 = flagO2, flagCH4 = flagCH4, flagN2O = flagN2O,
      KmO2_Rh = KmO2_Rh, alphaVmaxSx = alphaVmaxSx,   EaVmaxSx = EaVmaxSx,
      KMSx = KMSx, 
      alphaVmaxCH4prod = alphaVmaxCH4prod, 
      EaVmaxCH4prod = EaVmaxCH4prod,  KMSx_CH4prod = KMSx_CH4prod,
      alphaVmaxCH4ox = alphaVmaxCH4ox, EaVmaxCH4ox = EaVmaxCH4ox,
      KMCH4ox = KMCH4ox, KmO2_CH4ox = KmO2_CH4ox, kl_CH4prod = kl_CH4prod,
      
      alphaVmaxN2Oprod_nitrif = alphaVmaxN2Oprod_nitrif, 
      EaVmaxN2Oprod_nitrif = EaVmaxN2Oprod_nitrif,  KM_N2Oprod_nitrif = KM_N2Oprod_nitrif,
      KMO2_N2Oprod_nitrif = KMO2_N2Oprod_nitrif, alphaVmaxN2Oprod_denitrif = alphaVmaxN2Oprod_denitrif,
      EaVmaxN2Oprod_denitrif = EaVmaxN2Oprod_denitrif, KM_N2Oprod_denitrif = KM_N2Oprod_denitrif,
      alphaVmaxN2Ored = alphaVmaxN2Ored, EaVmaxN2Ored = EaVmaxN2Ored, KMN2Ored = KMN2Ored,
      kl_N2Oprod = kl_N2Oprod, kl_N2Ored = kl_N2Ored, KMSx_denitrif = KMSx_denitrif,
      
      reltol = reltol, verbose = verbose,
      sxtot.range = sxtot.range, nquantile = nquantile, mult=mult,
      soilm.cv = soilm.cv, soilm.range=soilm.range, soilm.nlevel = soilm.nlevel,
      lognorm = lognorm)
    
    tmp
  }
  
  ## storage 1:n is the n row observation
  for(i in 1:n){
    tmp <- big.lst[[i]]
    ## storage
    meanRH[i] <- tmp$value[1]
    meanCH4[i] <- tmp$value[2]
    meanN2O[i] <- tmp$value[3]
    flag[[i]] <- tmp$flag
    olst[[i]] <- cbind(Data.point=i,micro.site=1:nrow(tmp$output),tmp$microsite, tmp$output)
    O2.lst[[i]] <- tmp$output[,"O2"]
    CH4.lst[[i]] <- tmp$output[,"CH4"]
    N2O.lst[[i]] <- tmp$output[,"N2O"]
  }
  
  # ## marginalized log-likelihood
  # lik__ <- mlik(
  #   y=cbind(Rh__,CH4__),
  #   mu=cbind(meanRH,meanCH4),
  #   nu=nu,R0=R0)
  
  ## return
  list(data=cbind(SoilM__,SoilT__,Rh__,CH4__,N2O__,NH4__,NO3__),
       fit_Rh=meanRH,fit_CH4=meanCH4,fit_N2O=meanN2O,par=x,tpar=z,
       flag=do.call(cbind,flag),
       output=do.call(rbind,olst),
       O2.lst=O2.lst,
       CH4.lst=CH4.lst,
       N2O.lst=N2O.lst,
       elapsed=proc.time()-ts0)
}
#not run-----------------------------------
#need to add N2O---
test <- function()
{
  rm(list=ls())
  
  ## Session > Work Dir > Source File Loc
  trench <- read.csv("2015data.csv",head=T)
  #batch <- floor(0:(nrow(trench)-1)/10)+1
  #ii <- with(trench,which(batch<12 & batch %% 2==0))
  #trench$Rh[ii] <- NA
  #trench$CH4[ii] <- NA
  #write.csv(trench,file="mf13b_est.csv")
  
  ## Initial values to start estimation
  z0 <-   c(
    KmO2_Rh=245,alphaVmaxSx=230000000000, EaVmaxSx=72,
    KMSx=1000,
    alphaVmaxCH4prod=3000000,  EaVmaxCH4prod=100,
    KMSx_CH4prod=350,
    alphaVmaxCH4ox=7, EaVmaxCH4ox=30,
    KMCH4ox=1.00E-02,  KmO2_CH4ox=43,kl_CH4prod=3,
    alphaVmaxN2Oprod_nitrif=100, EaVmaxN2Oprod_nitrif=50,  KM_N2Oprod_nitrif=22.22,KMO2_N2Oprod_nitrif=143.3,
    alphaVmaxN2Oprod_denitrif=100, EaVmaxN2Oprod_denitrif=50,  KM_N2Oprod_denitrif=26,
    alphaVmaxN2Ored=1500, EaVmaxN2Ored=50,   KMN2Ored=0.16, 
    kl_N2Oprod=14.3,kl_N2Ored=7.2,KMSx_denitrif=1000
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
    kl_CH4prod=c(0.03,300,1),
    
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

  # source("damm_11a.R")
  source("damm_11a-non-steady-state.R")
  source("damm_utils_3.R")
 
  x0 <- damm_btran(z0,hyper=hyper)
  
  source("lulik_4b.R")
  debug(lulik)
  r0 <- lulik(
    x=x0,formula=~SoilM+SoilT+Rh+CH4+N2O+NH4+NO3,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209,  CH4airfrac = 1.8,N2Oairfrac=0.328,BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth_O2 = 10, Depth_CH4 = 10,Depth_N2O =10, microsite = 0,total.microsite = 10000, 
    time.step = 300,
    hyper=hyper
  )
  
  library(doParallel)
  cl <- makeCluster(detectCores()-1)
  registerDoParallel(cl)
  sink <- clusterEvalQ(cl,  source("damm_11a-non-steady-state.R"))
  sink <- clusterEvalQ(cl,source("damm_utils_3.R"))
  
  r1 <- lulik(
    x=x0,formula=~SoilM+SoilT+Rh+CH4+N2O+NH4+NO3,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209,  CH4airfrac = 1.8,N2Oairfrac=0.328,BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth_O2 = 10, Depth_CH4 = 10,Depth_N2O =10, microsite = 0,total.microsite = 10000, time.step = 300,
    hyper=hyper
  )
  
  stopCluster(cl)
  #plot(r0$fit_Rh,r1$fit_Rh,xlab="Serial",ylab="Parallel",main="Rh")
  #plot(r0$fit_CH4,r1$fit_CH4,xlab="Serial",ylab="Parallel",main="CH4")
  
  #r1$elapsed*length(cl)/r0$elapsed
  # user     system    elapsed 
  # 0.08473684 7.00000000 1.46963351 
  png(file="lulik_4b.png",width=3000,height=3000,res=400)
  op <- par(mar=c(4.1,4.1,1.1,1.1),mfrow=c(2,1))
  with(r1,{
    plot(data[,"Rh__"],bg="steelblue1",col=1,pch=21,cex=2,
         xlab="Time",ylab="Rh",main="",
         ylim=range(data[,"Rh__"],fit_Rh,na.rm=T))
    lines(fit_Rh,lty=2)
    plot(data[,"CH4__"],bg="salmon",col=1,pch=21,cex=2,
         xlab="Time",ylab="CH4",main="",
         ylim=range(data[,"CH4__"],fit_CH4,na.rm=T))
    lines(fit_CH4,lty=2)
  })
  par(op)
  dev.off()
}

