## optimizeP.R

## helper script to optimize parameters of the liquid steady state model 
## using observed time series
## of SoilM, SoilT, Rh and CH4 data and N2O data

## Arguments:
## z0 : initial parameter values to start MCMC
## hyper: initial parameter range and whether on log-scale for each parameter
## trench.file: excel data file name
## microsite.file: excel data file for fixed microsite
##                 0 to generate bivariate distribution automatically
## num.node      : number of cores for parallel computing,
##                default total number of cores minus one
## nsweep        : number of sweeps through the parameter space
## outfile       : root name for output files ending with
## tol           : relative tolerence for each conditional optimization
## ...           : additional argument to lulik
## value
##          : a time series plot of observed and fitted Rh and CH4 values
## param.csv: for final parameter values
## main.csv:  for data point and microsite specific output
## Rh.csv:   for data point specific Rh value
## CH4.csv:  for data point specific CH4 value
## N2O.csv:  for data point specific N2O value
#######################################################

## load required packages
require(truncnorm)
require(doParallel)

## Load R functions
source("damm_utils.R")
source("damm.R")
source("lulik.R")
source("lulik_utils.R")

optimizeP <- function(z0,hyper,trench.file,outfile,tol,microsite.file=0,num.node=NULL,nsweep=9,...)#nsweep=999
{
  ## log transformed initial parameters z0 in ranges (hyper) on the real line: z0--->x0 (hyper)
  x0 <- damm_btran(z=z0,hyper=hyper)
 
  ## Load data (observation data)
  trench <- read.csv(trench.file,head=T)
  
  ## prepare the high performance cluster
  if(is.null(num.node)){
    num.node <- detectCores() - 1 
  }
  cl <- makeCluster(num.node)
  registerDoParallel(cl)
  sink <- clusterEvalQ(cl,source("damm.R"))#damm, damm_work
  sink <- clusterEvalQ(cl,source("damm_utils.R"))#damm_btran damm_tran
  
  ## Storage for the algorithm
  history <- matrix(NA,nsweep,length(x0))#nsweep rows 25 column parameters
  fnval_CH4_N2O<- vector("list",nsweep)
  O2.val <- vector("list",nsweep)
  CH4.val <- vector("list",nsweep)
  N2O.val <- vector("list",nsweep)
  
  ## start the iterations with optimizing O2 while fixing CH4,N2O
  ####x0[1:4]: KmO2_Rh=245,alphaVmaxSx=230000000000,EaVmaxSx=72,KMSx=1000
  nls.O2.0 <- optim(
    par=x0[1:4],fn=lulik.rh.optim,gr=NULL,method="Nelder-Mead",control=list(reltol=tol),
    z=x0[5:25],formula=~SoilM+SoilT+Rh+CH4+N2O+NH4+NO3,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8, N2Oairfrac=0.328,
    BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10, Depth_N2O=10, microsite = 0,
    total.microsite = 10000, time.step = 300,
    iO2.lst=NULL,iCH4.lst=NULL,iN2O.lst=NULL,hyper=hyper#,...
  )
  ## evaluate loglik at the last values
  nls.O2.r0 <- lulik(
    x=c(nls.O2.0$par,x0[5:25]),formula=~SoilM+SoilT+Rh+CH4+N2O+NH4+NO3,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8,N2Oairfrac=0.328,
    BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10,Depth_N2O=10,microsite = 0,total.microsite = 10000, time.step = 300,
    flagO2=TRUE, flagCH4=FALSE, flagN2O=FALSE,hyper=hyper#,...
  )
  nls.O2.x0 <- c(nls.O2.0$par,x0[5:25])
  
  for(i in 1:nsweep){
    ## continue iteration with CH4,N2O  while fixing O2
    cat("iter ",i," CH4,N2O given O2.\n")
    #browser()
    nls.CH4.N2O <- optim(
      par=nls.O2.x0[5:25],fn=lulik.ch4.n2o.optim,gr=NULL,method="Nelder-Mead",control=list(reltol=tol),
      z=nls.O2.x0[1:4],formula=~SoilM+SoilT+Rh+CH4+N2O+NH4+NO3,data=trench,nlevels=100,
      R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8,N2Oairfrac=0.328,
      BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
      p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10,Depth_N2O=10,
      microsite = 0,total.microsite = 10000, time.step = 300,
      iO2.lst=nls.O2.r0$O2.lst,iCH4.lst=nls.O2.r0$CH4.lst,iN2O.lst=nls.O2.r0$N2O.lst,
      hyper=hyper#,...
    )
    ## evaluate loglik at the last values
    nls.CH4.N2O.r1 <- lulik(
      x=c(nls.O2.x0[1:4],nls.CH4.N2O$par),formula=~SoilM+SoilT+Rh+CH4+N2O+NH4+NO3,
      data=trench,nlevels=100,
      R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8,N2Oairfrac=0.328,
      BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
      p=0.00046,RQ=1, Depth_O2=10, Depth_CH4 = 10,Depth_N2O=10,
      microsite = 0,total.microsite = 10000, time.step = 300,
      iO2.lst=nls.O2.r0$O2.lst,iCH4.lst=nls.O2.r0$CH4.lst,iN2O.lst=nls.O2.r0$N2O.lst,
      hyper=hyper#,...
    )
    nls.CH4.N2O.x1 <- c(nls.O2.x0[1:4],nls.CH4.N2O$par)
    
    
    ## store results
    fnval_CH4_N2O[[i]] <- nls.CH4.N2O #The best set of parameters found.
    # fnval_N2O[[i]] <- nls.N2O
    
   # (history[i,] <- nls.CH4.r1$tpar) ## original scale parameters
    history[i,] <- nls.CH4.N2O.r1$tpar
    
    # O2.val[[i]] <- nls.CH4.r1$O2.lst
    # CH4.val[[i]] <- nls.CH4.r1$CH4.lst
    O2.val[[i]] <- nls.CH4.N2O.r1$O2.lst
    CH4.val[[i]] <- nls.CH4.N2O.r1$CH4.lst
    N2O.val[[i]] <- nls.CH4.N2O.r1$N2O.lst
    
    cat(history[i,],", convergence=",nls.CH4.N2O$convergence,", fn=",nls.CH4.N2O$value,"\n",
        file=paste(outfile,"_obj_CH4_N2O.txt",sep=""),append=TRUE)
    save(O2.val,CH4.val,N2O.val,history,i,file=paste(outfile,"_dump_CH4_N2O.rda",sep=""))
    
    ## repeat iteration with O2 while fixing CH4.N2O
    cat("iter ",i," O2 given CH4, N2O.\n")
    nls.O2.1 <- optim(
      par=nls.CH4.N2O.x1[1:4],fn=lulik.rh.optim,gr=NULL,method="Nelder-Mead",control=list(reltol=tol),
      z=nls.CH4.N2O.x1[5:25], formula=~SoilM+SoilT+Rh+CH4+N2O+NH4+NO3,
      data=trench,nlevels=100,
      R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8,N2Oairfrac=0.328,
      BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
      p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10,Depth_N2O=10,
      microsite = 0,total.microsite = 10000, time.step = 300,
      iO2.lst=nls.CH4.N2O.r1$O2.lst,iCH4.lst=nls.CH4.N2O.r1$CH4.lst,iN2O.lst=nls.CH4.N2O.r1$N2O.lst,
      hyper=hyper#,...
    )
    ## evaluate the latest likelihood function
    nls.O2.r1 <- lulik(
      x=c(nls.O2.1$par,nls.CH4.N2O.x1[5:25]),
      formula=~SoilM+SoilT+Rh+CH4+N2O+NH4+NO3,data=trench,nlevels=100,
      R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8,N2Oairfrac=0.328,
      BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
      p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10,Depth_N2O=10,
      microsite = 0,total.microsite = 10000, time.step = 300,
      flagO2=TRUE, flagCH4=FALSE,flagN2O=FALSE,
      iO2.lst=nls.CH4.N2O.r1$O2.lst,iCH4.lst=nls.CH4.N2O.r1$CH4.lst,iN2O.lst=nls.CH4.N2O.r1$N2O.lst,
      hyper=hyper#,...
    )
    nls.O2.x1 <- c(nls.O2.1$par,nls.CH4.N2O.x1[5:25])
    
    ## stop condition
    delta <- abs(nls.O2.x1-nls.O2.x0)
    cat(delta,"\n")
    if(all(delta < tol * abs(nls.O2.x0))){
      cat("Converged a iteration ",i,"\n")
      break
    }
    
    ## update
    nls.O2.x0 <- nls.O2.x1
  }
  
  ## prepare output
  x1 <- nls.O2.x1
  output <- lulik(
    x=x1, formula=~SoilM+SoilT+Rh+CH4+N2O+NH4+NO3,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209, CH4airfrac = 1.8,N2Oairfrac=0.328,
    BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10,Depth_N2O=10,
    microsite = 0,total.microsite = 10000, time.step = 300,
    hyper=hyper)#,...)
  save(fnval_CH4_N2O,history,O2.val,CH4.val,N2O.val,output, file=paste(outfile,".rda",sep=""))  ## save compressed R data
  write.csv(output$tpar, file=paste(outfile,"_parm.csv",sep=""),row.names = FALSE)
  write.csv(output$output,file=paste(outfile,"_main.csv",sep=""),row.names = FALSE)
  write.csv(output$fit_Rh, file=paste(outfile,"_Rh.csv",sep=""),row.names = FALSE)
  write.csv(output$fit_CH4, file=paste(outfile,"_CH4.csv",sep=""),row.names = FALSE)
  write.csv(output$fit_N2O, file=paste(outfile,"_N2O.csv",sep=""),row.names = FALSE)
  
  ## fit versus observed data
  png(file=paste(outfile,".png",sep=""),width=3000,height=3000,res=400)
  op <- par(mar=c(4.1,4.1,1.1,1.1),mfrow=c(3,1))
  with(output,{
    plot(data[,"Rh__"],bg="steelblue1",col=1,pch=21,cex=2,
         xlab="Time",ylab="Rh",main="",
         ylim=range(data[,"Rh__"],fit_Rh,na.rm=T))
    lines(fit_Rh,lty=2)
    plot(data[,"CH4__"],bg="salmon",col=1,pch=21,cex=2,
         xlab="Time",ylab="CH4",main="",
         ylim=range(data[,"CH4__"],fit_CH4,na.rm=T))
    lines(fit_CH4,lty=2)
    
    plot(data[,"N2O__"],bg="salmon",col=1,pch=21,cex=2,
         xlab="Time",ylab="N2O",main="",
         ylim=range(data[,"N2O__"],fit_N2O,na.rm=T))
    lines(fit_N2O,lty=2)
  })
  par(op)
  dev.off()
  ## stop cluster
  stopCluster(cl)
  
  output
}

test <- function()
{
  ## in model_framework_14_a.R
}
