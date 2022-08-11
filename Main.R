## model main function
rm(list=ls())
source("optimizeP.R")
source("simulate.R")

# model 3 input files: 
#"initial_par_input.csv","initial_parrange_input.csv","data.csv"
#read inital parameters from initial_par_input.csv
init_par <- read.csv("initial_par_input.csv",header = T)
z0 <- init_par[[2]]
names(z0) <- init_par[[1]]
#read initial range and whether on log-scale for each parameter
init_parrange <- read.csv("initial_parrange_input.csv",header = T,row.names = 1)
hyper <- as.list(init_parrange)

#note:there are 25 parameter initial and ranges for CH4 and N2O modules.

n.mcmc <- 99  ## number of MCMC approximate samples
##  re-sample the data using a gam model implemented in INLA for posterior
##  predictive analyses, initialized with warm starts
trench <- read.csv("data.csv",head=T) #data.csv 
trench$t__ <- 1:nrow(trench)
set.seed(1)
Rh.samples <- simulate(Rh~s(t__),data=trench,nsweep=n.mcmc)
set.seed(2)
N2O.samples <- simulate(N2O~s(t__),data=trench,nsweep=n.mcmc)
set.seed(3)
CH4.samples <- simulate(CH4~s(t__),data=trench,nsweep=n.mcmc)

## start the posterior predictive analyses
z1 <- z0
result <- vector("list",n.mcmc)
tmp <- trench
cat("Start Posterior Predictive Analyses\n")
for(i in 1:n.mcmc){
  cat("MC sample =",i,"/",n.mcmc,"\n")
  tmp$Rh <- Rh.samples[,i]
  tmp$CH4 <- CH4.samples[,i]
  tmp$N2O <- N2O.samples[,i]
  write.csv(tmp,file="data_tmp.csv")
  result[[i]] <- optimizeP(z0=z1,hyper=hyper,trench.file="data_tmp.csv",outfile="optim_parameter",
                              tol=5e-3)
}

## prepare output
parm <- sapply(result,function(elmt) elmt$tpar)
parm.out <- t(apply(parm,1,quantile,prob=c(0.5,0.025,0.975)))
names(parm.out) <- c("Estimate","LCL","UCL")
write.csv(parm.out,file="all_parm_cl.csv",row.names=TRUE)

Rh.fit <- sapply(result,function(elmt) elmt$fit_Rh)
N2O.fit <- sapply(result,function(elmt) elmt$fit_N2O)
CH4.fit <- sapply(result,function(elmt) elmt$fit_CH4)


## plot the fitted data.
png(file=paste("range_fit",".png",sep=""),width=3000,height=2500,res=400)
op <- par(mar=c(4.1,4.1,1.1,1.1),mfrow=c(3,1))
plot(ts(Rh.fit),plot.type="single",type="l",lwd=1.2,col="steelblue1",
     xlab="Time",ylab="Rh",main="",
     ylim=range(Rh.fit,trench$Rh,na.rm=T))
points(ts(trench$Rh),pch=21,cex=1.2,col=1,bg="steelblue1")
plot(ts(CH4.fit),plot.type="single",type="l",lwd=1.2,col="salmon",
     xlab="Time",ylab="CH4",main="",
     ylim=range(CH4.fit,trench$CH4,na.rm=T))
points(ts(trench$CH4),pch=21,cex=1.2,col=1,bg="salmon")
plot(ts(N2O.fit),plot.type="single",type="l",lwd=1.2,col="limegreen",
     xlab="Time",ylab="N2O",main="",
     ylim=range(N2O.fit,trench$N2O,na.rm=T))
points(ts(trench$N2O),pch=21,cex=1.2,col=1,bg="limegreen")

par(op)
dev.off()
