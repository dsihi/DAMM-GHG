
## ver. 7: allow formula and data to be used directly for INLA
## bug: no smooth interface.
## helper function to define smooth term fitting in INLA
## supporting offset through formula instead
require(mgcv)
source("getVars.R")
## helper function to fix smooth terms in INLA
## formula: the formula containing only smooth terms
## data: the data to define smooth terms
## U  : standard deviations of T Normal priors on standard deviation of random effects
##    : a vector of two, the first entry is for the smooth with rank > 1
##    :                  the second entry is for the smooth with rank 1
##    : NOTE: these numbers appear work with single smooth, but the linear
##    combination part is tricky with z models
## rankone: logical, whether include rankone smooth or not.
## family : glm family to process
## zmat   : character string for zmatrix list name
## value: X: the Z matrices for INLA definition
##        formula: the INLA formula including smooth terms
##        data: the expanded data frame for INLA model
##        sterms: smooth terms
##        zmat  : input for zmat name
s2inla <- function(formula,data,U=c(0.00001,10),rankone=FALSE,family=gaussian,zmat="X")
{
  l.tt <- terms(formula)
  l.ss <- grep("^s(.+)",attr(l.tt,"term.labels"))
  l.ff <- grep("^f(.+)",attr(l.tt,"term.labels"))
  l.oo <- row.names(attr(l.tt,"factors"))[attr(l.tt,"offset")]
  l.terms.nos <- attr(l.tt,"term.labels")[-l.ss] ## identify non-smooth terms
  l.response <- getVars(formula)[1]              ## identify response variable
  
  if(length(l.ss)==0){
    ## the formula contain no s term
    catch <- try(model.matrix(formula,data),silent = TRUE)
    frame <- model.frame(formula,data)
    offset <- model.offset(frame)
    if(is.null(offset)){
      offset <- rep(0,length(model.response(frame)))
    }
    return(list(y=model.response(frame),offset=offset,
                Design=catch,data=data,formula=formula,sterms=NULL,zmat=zmat))
  }
  
  ## remove the f portion of the formula to be compatible with jagam
  l.formula <- formula
  if(length(l.ff)>0){
    l.formula <- reformulate(c(getVars(formula)[-1][-c(l.ff,l.ss)],
                               attr(l.tt,"term.labels")[l.ss]),
                             response=l.response)
  }
  #browser()
  
  ## remove missing values
  data <- na.omit(data[,all.vars(formula)])
  jf <- tempfile()
  jd <- jagam(l.formula,data=data,file=jf,
              diagonalize = TRUE,na.action=na.pass,family=family)
  
  ## creat fixed design matrix
  nms <- names(jd$pregam$cmX) ## if this name is not null, it is a design matrix
  Design <- jd$pregam$X[,nchar(nms)>0] ## now without f terms
  #browser()
  offset <- model.offset(model.frame(jd$pregam$terms,data))
  if(is.null(offset)){
    offset <- rep(0,nrow(data))
  }
  
  ## create list for random design matrix, formula and simply IDs
  X <- vector("list",length(jd$pregam$smooth))
  fout <- vector("list",length(X))
  varnames <- letters[1:length(X)]
  ones <- vector("list",length(X))
  
  for(i in 1:length(jd$pregam$smooth)){ ## loop thru smooth term
    ## create id matrices
    tmp <- matrix(1:nrow(data),nrow(data),length(jd$pregam$smooth[[i]]$rank))
    colnames(tmp) <- paste(varnames[i],1:length(jd$pregam$smooth[[i]]$rank),sep="")
    ones[[i]] <- as.data.frame(tmp)
    
    ## create design matrices
    X[[i]] <- vector("list",length(jd$pregam$smooth[[i]]$rank))
    start <- jd$pregam$smooth[[i]]$first.para
    ff <- rep("",length(jd$pregam$smooth[[i]]$rank))
    for(j in 1:length(jd$pregam$smooth[[i]]$rank)){
      idx <- start+seq(0,jd$pregam$smooth[[i]]$rank[j]-1)
      cat("idx=",idx,"\n")
      if(length(idx)>1){ ## smooth rank > 1
        l.U <- U[1]
      }                  ## end
      else{              ## smooth rank == 1
        l.U <- U[2]
      }                  ## end
      X[[i]][[j]] <- jd$jags.data$X[,idx,drop=FALSE]
      ## INLA formula
      tmpfor <- paste("f(",varnames[i],j,',model="z",',"Z=",zmat,"[[",i,"]][[",j,"]],",
                      'hyper=list(prec=list(prior="logtnormal",param=c(0,',
                      l.U,'))))',
                      sep="")
      if(length(idx)>1){
        ff[j] <- tmpfor
      }
      else{
        if(rankone){
          ff[j] <- tmpfor
        }
      }
      start <- max(idx)+1
    }
    fout[[i]] <- paste(ff[nchar(ff)>0],collapse = "+")
  }
  
  #browser()
  ## process data for INLA
  r.data <- cbind(data,do.call(cbind,ones))
  ## process formula for INLA
  
  r.formula <- reformulate(c(do.call(c,fout),l.terms.nos,l.oo),response=l.response)
  l.r <- list(y=jd$pregam$y,Design=Design,X=X,formula=r.formula,data=r.data,
       pregam=jd$pregam,offset=offset,
       sterms=getVars(formula)[-1][l.ss],zmat=zmat)
  names(l.r)[3] <- zmat
  
  l.r
}

test7 <- function(){
  rm(list=ls())
  gc(T)
  #setwd("C:/Users/dliang/Google Drive/Summer_2017/Turtle_Watch/Analysis")
  require(mgcv)
  set.seed(2) ## simulate some data... 
  n <- 400
  dat <- gamSim(1,n=n,dist="normal",scale=2)
  dat$y[1] <- NA
  dat$x0[2] <- NA
  
  source("s2inla_7.R")
  #debug(s2inla)
  #sout <- s2inla(y~x0+f(inla.group(x1),model='rw2')+s(x2)+offset(log(x3)),data=dat)
  #sout <- s2inla(y~x0+s(x1)+s(x2)+offset(log(x3)),data=dat)
  sout <- s2inla(y~x0+s(x1)+s(x2),data=dat)
  # ## for smooth terms, the smooth and rank deficieny terms size
  # ## the rank deficiency term was not included in modeling.
  # ## select non-rank deficient terms
  # l.n <- as.vector(sapply(sout$X,function(elmt){
  #   sapply(elmt,ncol)
  # }))
  # l.id <- rep(1:length(l.n),l.n)
  # l.sel <- which(l.id %% 2==1)

  #X <- sout$X
  #for1 <- formula(paste("y~x2+x3+",paste(sout$formula,collapse="+"),sep=""))
  library(INLA)
  assign("X",sout$X)
  result <- inla(sout$formula,data=sout$data,
                 control.inla=list(int.strategy='auto'))
  
  b0 <- gam(y~x0+s(x1)+s(x2),data=dat,method="REML")
  plot(b0,select=1)
  n1 <- nrow(sout$data)
  points(sout$data$x1,result$summary.random$a1[1:n1,5])
  points(sout$data$x1,result$summary.random$a1[1:n1,4],cex=0.5)  
  points(sout$data$x1,result$summary.random$a1[1:n1,6],cex=0.5)
  
  plot(b0,select=2)
  points(sout$data$x2,result$summary.random$b1[1:n1,5])
  points(sout$data$x2,result$summary.random$b1[1:n1,4],cex=0.5)  
  points(sout$data$x2,result$summary.random$b1[1:n1,6],cex=0.5)
  
  summary(result)$fixed
  summary(b0)
}

test8 <- function(){
  rm(list=ls())
  gc(T)
  setwd("C:/Users/dliang/Google Drive/Summer_2017/Turtle_Watch/Analysis")
  require(mgcv)
  set.seed(2) ## simulate some data... 
  n <- 400
  dat <- gamSim(1,n=n,dist="normal",scale=2)
  dat$y[1] <- NA
  dat$x0[2] <- NA
  
  source("s2inla_7.R")
  #debug(s2inla)
  #sout <- s2inla(y~x0+f(inla.group(x1),model='rw2')+s(x2)+offset(log(x3)),data=dat)
  #sout <- s2inla(y~x0+s(x1)+s(x2)+offset(log(x3)),data=dat)
  sout <- s2inla(y~x0+s(x1)+s(x2)+offset(log(x3)),data=dat)
  # ## for smooth terms, the smooth and rank deficieny terms size
  # ## the rank deficiency term was not included in modeling.
  # ## select non-rank deficient terms
  # l.n <- as.vector(sapply(sout$X,function(elmt){
  #   sapply(elmt,ncol)
  # }))
  # l.id <- rep(1:length(l.n),l.n)
  # l.sel <- which(l.id %% 2==1)
  
  #X <- sout$X
  #for1 <- formula(paste("y~x2+x3+",paste(sout$formula,collapse="+"),sep=""))
  library(INLA)
  assign(sout$zmat,sout$X)
  result <- inla(sout$formula,data=sout$data,
                 control.inla=list(int.strategy='auto'))
  result$summary.smooth <- result$summary.random
  names(result$summary.smooth) <- sout$sterms
  str(result$summary.smooth)
  
  b0 <- gam(y~x0+s(x1)+s(x2)+offset(log(x3)),data=dat,method="REML")
  plot(b0,select=1)
  n1 <- nrow(sout$data)
  points(sout$data$x1,result$summary.random$a1[1:n1,5])
  points(sout$data$x1,result$summary.random$a1[1:n1,4],cex=0.5)  
  points(sout$data$x1,result$summary.random$a1[1:n1,6],cex=0.5)
  
  plot(b0,select=2)
  points(sout$data$x2,result$summary.random$b1[1:n1,5])
  points(sout$data$x2,result$summary.random$b1[1:n1,4],cex=0.5)  
  points(sout$data$x2,result$summary.random$b1[1:n1,6],cex=0.5)
  
  summary(result)$fixed
  summary(b0)
}

