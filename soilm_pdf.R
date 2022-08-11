## soilm_pdf_3.R
## build microsite bivariate distribution

## helper script to build sxtot distribution given microsite
# sxtot.range  assumed lower and upper limit for Co2
# nquantile    number of levels for sxtot
# nsite        total number of microsites
# mult         multiplier between 0 and 1
##             the smaller the more skewed the output defaut 1/4
## value: sxtot and freqency
sxtot.microsite <- function(sxtot.range,nquantile,nsite,mult=1/4,verbose=TRUE)
{
  ## back-calculate the log-normal distribution using 2 sigma rule
  sxtot.log.mean <- mean(log(sxtot.range))
  sxtot.log.sd  <- diff(log(sxtot.range)) * mult
  
  ## equally divide the range according to the number of levels
  sxtot.q <- seq(sxtot.range[1],sxtot.range[2],len=nquantile-1)
  sxtot.mid <- (sxtot.q[-1]+sxtot.q[-length(sxtot.q)])/2
  sxtot.prob.vec <- pnorm(log(sxtot.q),mean=sxtot.log.mean,sd=sxtot.log.sd)
  sxtot.prob <- diff(sxtot.prob.vec)
  ## and calculate the proportions
  
  ## return the microsite distribution
  sxtot.microsite <- data.frame(
    sxtot=sort(c(sxtot.range,sxtot.mid)),
    freq=nsite*c(pnorm(log(sxtot.range[1]),sxtot.log.mean,sxtot.log.sd),sxtot.prob,
                 pnorm(log(sxtot.range[2]),sxtot.log.mean,sxtot.log.sd,lower.tail=FALSE)))
  
  sxtot.microsite
}

test1 <- function()
{
  rm(list=ls())
  ## Session > Work Dir > Source File Loc
  source("soilm_pdf_3.R")
  plot(sxtot.microsite(c(0.01,0.15),5,1000,mult=1/4),ylim=c(0,1000))
  plot(sxtot.microsite(c(0.01,0.15),5,1000,mult=1/5),ylim=c(0,1000))
  plot(sxtot.microsite(c(0.01,0.15),5,1000,mult=1/10),ylim=c(0,1000))
}

## helper script to build soil moisture distribution
## soilm          observed soil moisture to calibrate
## soilm.cv       assumed coefficient of variation (%) for soil moisture
## soilm.range    assumed lower and upper range (%) for soil moisture
## soilm.nlevel   number of levels to develop soil moisture distribution
## lognorm        whether lognormal or normal distribution
require(truncnorm)
## value: soilm and probability
soilm.microsite <- function(soilm,soilm.cv,soilm.range,soilm.nlevel,lognorm=TRUE){
  if(lognorm){
    soilm.limit <- log(soilm.range*soilm/100)
    soilm <- log(soilm)
    soilm.sd <- sqrt(log(1+(soilm.cv/100)^2))
    soilm.q <- seq(soilm.limit[1],soilm.limit[2],length=soilm.nlevel+1)
    soilm.mid <- (soilm.q[-1]+soilm.q[-length(soilm.q)])/2
    soilm.prob.vec <- ptruncnorm(soilm.q,a=soilm.limit[1],b=soilm.limit[2],mean=soilm,sd=soilm.sd)
    soilm.prob <- diff(soilm.prob.vec)
    soilm.microsite <- data.frame(
      soilm=exp(soilm.mid),prob=soilm.prob
    )
  }
  else{
    ## normal distribution
    soilm.sd <- soilm * soilm.cv/100
    soilm.q <- seq(soilm.range[1],soilm.range[2],length=soilm.nlevel+1)*soilm/100
    soilm.mid <- (soilm.q[-1]+soilm.q[-length(soilm.q)])/2
    soilm.prob.vec <- ptruncnorm(soilm.q,a=soilm.range[1]*soilm/100,
                                 b=soilm.range[2]*soilm/100,mean=soilm,sd=soilm.sd)
    soilm.prob <- diff(soilm.prob.vec)
    soilm.microsite <- data.frame(
      soilm=soilm.mid,prob=soilm.prob
    )
  }
  soilm.microsite
}

test2 <- function()
{
  rm(list=ls())
  ## Session > Work Dir > Source File Loc
  source("soilm_pdf_3.R")
  plot(soilm.microsite(20,20,c(70,200),20),type="l")
}

## script to build bivariate sxtot and soilm distribution
soilm.pdf <- function(sxtot.microsite,soilm.microsite,file=NULL){
  distribution <- expand.grid(soilm.id=1:nrow(soilm.microsite),
                              sxtot.id=1:nrow(sxtot.microsite))
  distribution$sxtot <- sxtot.microsite$sxtot[distribution$sxtot.id]
  distribution$soilm <- soilm.microsite$soilm[distribution$soilm.id]
  distribution$freq <- sxtot.microsite$freq[distribution$sxtot.id] * soilm.microsite$prob[distribution$soilm.id]
  if(!is.null(file)){
    write.csv(distribution[,c("sxtot","soilm","freq")],file=file,row.names=FALSE)
  }
  distribution
}

test3 <- function()
{
  rm(list=ls())
  ## Session > Work Dir > Source File Loc
  source("soilm_pdf_3.R")
  soilm1 <- soilm.microsite(20,20,c(70,200),5)
  sxtot1 <- sxtot.microsite(c(0.01,0.15),5,10000,mult=1/10)
  dist1 <- soilm.pdf(sxtot1,soilm1,"microsite_5.csv")
}