#Non-steady state model for CH4-Liq

source("soilm_pdf.R") ## load script for bi-variate distribution

## Arguments
damm <- function( 
  R, O2airfrac, CH4airfrac, N2Oairfrac,BD, PD,  Porosity,
  p, RQ, 
  
  SoilM.o,   SoilT ,  NH4, NO3, ###added NO3 and NH4
  Depth_O2, Depth_CH4, Depth_N2O,#Depth ,   ###correct for diff depth
  #SoilM.o is soil moisture observation
  microsite,total.microsite,
  
  time.step,
  
  iO2, iCH4, iN2O,## initial O2 and CH4 values
  
  flagO2,  flagCH4, flagN2O,## whether run O2 and CH4 to steady states
  
  KmO2_Rh, alphaVmaxSx,   EaVmaxSx,   KMSx,   
  alphaVmaxCH4prod, EaVmaxCH4prod,  KMSx_CH4prod,
  alphaVmaxCH4ox, EaVmaxCH4ox,   KMCH4ox, KmO2_CH4ox, kl_CH4prod,
  alphaVmaxN2Oprod_nitrif, EaVmaxN2Oprod_nitrif,  KM_N2Oprod_nitrif, KMO2_N2Oprod_nitrif,
  alphaVmaxN2Oprod_denitrif, EaVmaxN2Oprod_denitrif,  KM_N2Oprod_denitrif,
  alphaVmaxN2Ored, EaVmaxN2Ored, KMN2Ored, 
  kl_N2Oprod, kl_N2Ored,KMSx_denitrif,
  
  reltol,verbose,
  sxtot.range=NULL, nquantile=NULL, mult=NULL,
  soilm.cv=NULL, soilm.range=NULL, soilm.nlevel=NULL, lognorm=NULL
)
{
  #First creat microsite at each obervation
  ## Build microsite bi-variate distribution 
  if(microsite==0){
    ## Check site specific limit
    if(SoilM.o*soilm.range[2]/100 > (Porosity*100)){
      soilm.range[2] <- (Porosity*100)*95/SoilM.o
      if(verbose>3){
        cat("Soil Moisture upper limit updated to ",soilm.range[2],"\n")
      }
    }

    soilm1 <- soilm.microsite(
      soilm=SoilM.o,soilm.cv=soilm.cv,soilm.range = soilm.range,
      soilm.nlevel = soilm.nlevel,lognorm=lognorm)
    sxtot1 <- sxtot.microsite(
      sxtot.range = sxtot.range,nquantile = nquantile,nsite=total.microsite,
      mult=mult,verbose = FALSE)
    microsite <- soilm.pdf(sxtot1,soilm1)
  }
  #
  #create microsite combination of soil carbon and soil moisture
  Sxtot <- microsite$sxtot
  SoilM <- microsite$soilm
  Freq <- microsite$freq

  ## initial values
  if(length(iO2)==1){
    iO2 <- rep(iO2,length(Sxtot))
  }
  if(length(iCH4)==1){
    iCH4 <- rep(iCH4,length(Sxtot))
  }
  if(length(iN2O)==1){
    iN2O <- rep(iN2O,length(Sxtot))
  }

  #### Need to check for possible changes to be done

  ## per micro-site steady state outputs and interim results
  vRh <- rep(NA,length(Sxtot))
  O2 <- rep(NA,length(Sxtot))
  vCH4 <- rep(NA,length(Sxtot))
  vN2O <- rep(NA,length(Sxtot))
  Rh_calc <- rep(NA,length(Sxtot))
  O2_calc <- rep(NA,length(Sxtot))
  CH4_diff <- rep(NA,length(Sxtot))
  CH4_calc <- rep(NA,length(Sxtot))
  N2O_diff <- rep(NA,length(Sxtot))
  N2O_calc <- rep(NA,length(Sxtot))
  a_calc <- rep(NA,length(Sxtot)) #need to remove
  flag <- rep(NA,length(Sxtot)) #need to remove
  output <- matrix(NA,length(Sxtot),28) ###CORRECT HERE to 18 # 16 without D_P and flag

  #Need to change here
  colnames(output) <- c("O2","O2_l_micromolperliter","O2_l_micromolpermicrosite","Rh_l_micromolpermicrosite","Rh",
                        "O2_l_micromolpermicrosite_conc","O2_l_micromolperliter_conc","CH4","CH4_l_micromolperliter",
                        "CH4_l_micromolpermicrosite","CH4prod","CH4_l_conc_micromolpermicrosite",
                        "CH4_l_conc_micromolperliter","CH4ox","CH4Flux_l_micromolpermicrosite",
                     "CH4Flux_l_micromolperliter","D_P","flag",
                     "N2O","N2O_l_micromolperliter",
                     "N2O_l_micromolpermicrosite","N2Oprod_nitrif","N2Oprod_denitrif","N2O_l_conc_micromolpermicrosite",
                     "N2O_l_conc_micromolperliter","N2Ored","N2OFlux_l_micromolpermicrosite",
                     "N2OFlux_l_micromolperliter")
    #loop each microsites (25 for each observation,different soil carbon and soil moisture)                  
    for(i in 1:length(Sxtot)){
    damm_return <- damm_work(
      R=R, O2airfrac=O2airfrac,CH4airfrac=CH4airfrac,N2Oairfrac=N2Oairfrac,BD=BD,PD=PD,
      Porosity=Porosity,p=p,RQ=RQ,
      SoilM.o=SoilM.o,SoilT=SoilT,Depth_O2=Depth_O2,Depth_CH4=Depth_CH4,Depth_N2O=Depth_N2O,#Depth=Depth, ###correct for diff depth
      Sxtot=Sxtot[i], SoilM=SoilM[i],
      total.microsite=total.microsite,NO3=NO3,NH4=NH4,
      time.step=time.step,
      iO2=iO2[i], iCH4=iCH4[i], iN2O=iN2O[i], ## initial O2 and CH4 values
      flagO2=flagO2,  flagCH4=flagCH4, flagN2O=flagN2O, ## whether run O2 and CH4 to steady states
      KmO2_Rh=KmO2_Rh, alphaVmaxSx=alphaVmaxSx,   EaVmaxSx=EaVmaxSx,
      KMSx=KMSx,alphaVmaxCH4prod=alphaVmaxCH4prod,
      EaVmaxCH4prod=EaVmaxCH4prod,KMSx_CH4prod=KMSx_CH4prod,
      alphaVmaxCH4ox=alphaVmaxCH4ox, EaVmaxCH4ox=EaVmaxCH4ox,
      KMCH4ox=KMCH4ox, KmO2_CH4ox=KmO2_CH4ox, kl_CH4prod=kl_CH4prod,
      
      alphaVmaxN2Oprod_nitrif=alphaVmaxN2Oprod_nitrif,EaVmaxN2Oprod_nitrif=EaVmaxN2Oprod_nitrif,
      KM_N2Oprod_nitrif=KM_N2Oprod_nitrif,KMO2_N2Oprod_nitrif=KMO2_N2Oprod_nitrif,
      alphaVmaxN2Oprod_denitrif=alphaVmaxN2Oprod_denitrif, EaVmaxN2Oprod_denitrif=EaVmaxN2Oprod_denitrif,  
      KM_N2Oprod_denitrif=KM_N2Oprod_denitrif,
      alphaVmaxN2Ored=alphaVmaxN2Ored, EaVmaxN2Ored=EaVmaxN2Ored,KMN2Ored=KMN2Ored, 
      kl_N2Oprod=kl_N2Oprod, kl_N2Ored=kl_N2Ored,KMSx_denitrif=KMSx_denitrif,
      
      reltol=reltol,verbose=verbose)
    output[i,] <- damm_return
    vRh[i] <- damm_return[5] #5
    O2[i] <- damm_return[1] #1
    vCH4[i] <- damm_return[16] #16 for CH4Flux_l_micromolperliter
    vN2O[i] <- damm_return[28] #28 for N2OFlux_l_micromolperliter
    
    a_calc[i] <- damm_return[17] #17
    flag[i] <- damm_return[18] # 18
    Rh_calc[i] <- vRh[i]*Freq[i]
    CH4_calc[i] <- vCH4[i] * Freq[i]
    N2O_calc[i] <- vN2O[i] * Freq[i]
    if(verbose>0) {
      cat(" micro site ",i," Rh=",vRh[i],"O2=",O2[i],"CH4=",vCH4[i]," flag=",flag[i],"\n\n")
    }

  }
  microsite.vol <- (3.14*12.5*12.5*10)/(total.microsite)
  Henry_constant_CO2 <- 3.4*10^(-2)*exp(2400*((1/(273.15+SoilT))-1/298.15)) ## derived constant and parameters
  Henry_constant_CH4 <- 1.4*10^(-3)*exp(1700*((1/(273.15+SoilT))-1/298.15)) ## derived constant and parameters
  Henry_constant_N2O <- 2.5*10^(-2)*exp(2600*((1/(273.15+SoilT))-1/298.15)) ## derived constant and parameters
  
  damm_final_return <- rep(NA,3)

  damm_final_return[1] <- (sum(Rh_calc)/sum(Freq))*(1/Henry_constant_CO2)*(10^-6)
  damm_final_return[1] <- damm_final_return[1]*12*(1/1000)*1000*10*10000*(3600/time.step)

  damm_final_return[2] <- (sum(CH4_calc)/sum(Freq))*(1/Henry_constant_CH4)*(10^-6)
  damm_final_return[2] <- damm_final_return[2]*12*(1/1000)*1000*10*10000*(3600/time.step)
  
  damm_final_return[3] <- (sum(N2O_calc)/sum(Freq))*(1/Henry_constant_N2O)
  damm_final_return[3] <- damm_final_return[3]*28*(1/1000)*10*10000*(3600/time.step)
  
  
  list(value=damm_final_return,flag=flag,output=output,microsite = microsite)
}

test2 <- function()
{
  rm(list=ls())
  ## Session > Set Work Dir > Source File Loc
  source("damm_11a-non-steady-state.R")
 # microsite <- read.csv("microsite.csv")
  microsite <- read.csv("microsite_5.csv")
  microsite <- 0
  sink("./damm_liq_debug_1.txt")
  #debug(damm_work)
  r <- damm( 
    R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8, N2Oairfrac=0.328,
    BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6, 
    p=0.00046,RQ=1,
    
    SoilM.o=30,   SoilT=20,   
    Depth_O2=10,Depth_CH4=10,Depth_N2O=10,#Depth=10 ,  
    microsite=microsite,  total.microsite = 10000, 
    NO3=0.7, NH4=5, #added NO3 and NH4
    time.step = 300,
    
    # iO2=rep(0.209,nrow(microsite)), 
    # iCH4=rep(1.8,nrow(microsite)),
    # iN2O=rep(0.328,nrow(microsite)),
    iO2=rep(0.209,1), 
    iCH4=rep(1.8,1),
    iN2O=rep(0.328,1),
  
    flagO2=TRUE, flagCH4=TRUE,flagN2O=TRUE,
    
    KmO2_Rh=245, alphaVmaxSx=230000000000,   EaVmaxSx=72,   KMSx=1000,   
    alphaVmaxCH4prod=3000000, EaVmaxCH4prod=100,  KMSx_CH4prod=350,
    alphaVmaxCH4ox=0.07, EaVmaxCH4ox=30,   KMCH4ox=1.00E-02, KmO2_CH4ox=43, kl_CH4prod=3,
    
    alphaVmaxN2Oprod_nitrif=100, EaVmaxN2Oprod_nitrif=50,  KM_N2Oprod_nitrif=22.22,KMO2_N2Oprod_nitrif=143.3,
    alphaVmaxN2Oprod_denitrif=100, EaVmaxN2Oprod_denitrif=50,  KM_N2Oprod_denitrif=26,
    alphaVmaxN2Ored=1500, EaVmaxN2Ored=50,   KMN2Ored=0.16, 
    kl_N2Oprod=14.3,kl_N2Ored=7.2,KMSx_denitrif=1000,
    reltol=sqrt(.Machine$double.eps),verbose=99,
    
    sxtot.range=c(0.01,0.15), nquantile=10, mult=1/4,
    soilm.cv=20, soilm.range=c(70,200), soilm.nlevel=10,lognorm=TRUE
  )
  sink()
  r$value
  #[1] 340.74690813  -0.05487782
  # excel
  #324, -0.05372
}

## damm_work:
## individual steady state for the given bivariate distribution of 
## Sxtot and soil moisture
## New Arguments added from version 9
## BD, PD
## KmO2_Rh (need optimization)
## Vmax_CH4prod (constant-need optimization?)
## KMSx_CH4prod (need optimization)
## Vmax_CH4ox (constant- need optimization?)
## KmO2_CH4ox (changed location from constant to parameter)
## kl_CH4prod (need optimization)
damm_work <- function( 
  R, O2airfrac,  CH4airfrac, N2Oairfrac,BD, PD,  Porosity,
  p, RQ, 
  
  SoilM.o,   SoilT ,   Sxtot, SoilM, Depth_O2,Depth_CH4,Depth_N2O,#Depth ,  ###correct for diff depth
  
  total.microsite,NO3, NH4,
  
  time.step,
  
  iO2, iCH4, iN2O,## initial O2 and CH4 values
  
  flagO2,  flagCH4, flagN2O,## whether run O2 and CH4 to steady states
  
  KmO2_Rh, alphaVmaxSx,   EaVmaxSx,   KMSx,   
  alphaVmaxCH4prod, EaVmaxCH4prod,  KMSx_CH4prod,
  alphaVmaxCH4ox, EaVmaxCH4ox,   KMCH4ox, KmO2_CH4ox, kl_CH4prod,
  
  alphaVmaxN2Oprod_nitrif, EaVmaxN2Oprod_nitrif,  KM_N2Oprod_nitrif,KMO2_N2Oprod_nitrif, 
  alphaVmaxN2Oprod_denitrif, EaVmaxN2Oprod_denitrif,  KM_N2Oprod_denitrif,
  alphaVmaxN2Ored, EaVmaxN2Ored, KMN2Ored, 
  kl_N2Oprod,kl_N2Ored,KMSx_denitrif,
  
  reltol,verbose
)
{
  if(verbose>3){
    cat("Soil M observed=",SoilM.o,", Soil T=",SoilT,
        "Sxtot =",Sxtot, ",Soil M =",SoilM,"\n")
  }
  if(verbose>3) {
    cat("KmO2_Rh=",KmO2_Rh," ",
        "alphaVmaxSx=",alphaVmaxSx, " ", 
        "EaVmaxSx=",EaVmaxSx, " ", 
        "KMSx=",KMSx, " ",
        "alphaVmaxCH4prod=",alphaVmaxCH4prod," ",
        "EaVmaxCH4prod=", EaVmaxCH4prod, " ",
        "KMSx_CH4prod", KMSx_CH4prod, " ",
        "alphaVmaxCH4ox=",alphaVmaxCH4ox," ", 
        "EaVmaxCH4ox=",EaVmaxCH4ox," ",
        "KMCH4ox=",KMCH4ox," ",
        "KmO2_CH4ox", KmO2_CH4ox, " ",
        "kl_CH4prod", kl_CH4prod, " ", ###Include Depth_O2 and Depth_CH4
        
        "alphaVmaxN2Oprod_nitrif=",alphaVmaxN2Oprod_nitrif," ",
        "EaVmaxN2Oprod_nitrif=", EaVmaxN2Oprod_nitrif, " ",
        "KM_N2Oprod_nitrif", KM_N2Oprod_nitrif, " ",
        "KMO2_N2Oprod_nitrif, ", KMO2_N2Oprod_nitrif, " ",
        "alphaVmaxN2Oprod_denitrif=",alphaVmaxN2Oprod_denitrif," ",
        "EaVmaxN2Oprod_denitrif=", EaVmaxN2Oprod_denitrif, " ",
        "KM_N2Oprod_denitrif", KM_N2Oprod_denitrif, " ",
        "alphaVmaxN2Ored=",alphaVmaxN2Ored," ", 
        "EaVmaxN2Ored=",EaVmaxN2Ored," ",
        "KMN2Ored=",KMN2Ored," ",
        "kl_N2Oprod", kl_N2Oprod, " ",
        "kl_N2Ored", kl_N2Ored, " ",
        "KMSx_denitrif=",KMSx_denitrif, " ", ###Include Depth_O2 and Depth_N2O
        "\n") 
  }
  
  #browser()
  ## derived constants 
  Dliq <- 1/Porosity^3 
  ## derived parameters
  VmaxCH4prod <- (alphaVmaxCH4prod)*exp(-EaVmaxCH4prod/(R*(SoilT+273.15))) 
  VmaxCH4ox <- (alphaVmaxCH4ox)*exp(-EaVmaxCH4ox/(R*(SoilT+273.15))) 
  
  VmaxN2Oprod_nitrif <- (alphaVmaxN2Oprod_nitrif)*exp(-EaVmaxN2Oprod_nitrif/(R*(SoilT+273.15))) 
  VmaxN2Oprod_denitrif <- (alphaVmaxN2Oprod_denitrif)*exp(-EaVmaxN2Oprod_denitrif/(R*(SoilT+273.15)))
  VmaxN2Ored <- (alphaVmaxN2Ored)*exp(-EaVmaxN2Ored/(R*(SoilT+273.15))) 
  
  microsite.vol <- (3.14*12.5*12.5*10)/(total.microsite)
  microsite.area <- (3.14*12.5*12.5)/(total.microsite)
  
  Henry_constant_O2 <- 1.3*10^(-3)*exp(1700*((1/(273.15+SoilT))-1/298.15)) ## derived constant and parameters
  Henry_constant_CO2 <- 3.4*10^(-2)*exp(2400*((1/(273.15+SoilT))-1/298.15)) ## derived constant and parameters
  Henry_constant_CH4 <- 1.4*10^(-3)*exp(1700*((1/(273.15+SoilT))-1/298.15)) ## derived constant and parameters
  Henry_constant_N2O <- 2.5*10^(-2)*exp(2600*((1/(273.15+SoilT))-1/298.15)) ## derived constant and parameters
  
  ## derived constant and parameters
  #D_P <-(2*0.261^3+0.04*0.261)*((Porosity-SoilM/100)/0.261)^(2+(3/2.9))
  D_P <- ((Porosity-SoilM/100)^(4.0/3.0))**(((SoilT+273.15)/293.15)^1.75)##a4_3 #update temp sensitivity in steadystate model also
  a4_3_max <- ((Porosity-0/100)^(4/3))*(((20+273.15)/293.15)^1.75)
  Dgas <- 1/a4_3_max
  DCH4 <- Dgas
  DO2 <- Dgas 
  DN2O <- Dgas #put correct number (and it's temp sensitivity) at STP from the paper Eric sent in steadystate model

  diffused_NO3 <- NO3*Dliq*(SoilM/100)^(3.0) ## derived constant and parameters
  diffused_NH4 <- NH4*Dliq*(SoilM/100)^(3.0) ## derived constant and parameters
  
  kl_CH4prod <- kl_CH4prod
  kMO2_Rh <- KmO2_Rh
  KMCH4ox <- KMCH4ox
  kMO2_CH4ox <- KmO2_CH4ox
  
  Sxtot_micromolperliter <- Sxtot*p*(100/SoilM)*1000*(10^6)*(1/12)
  Sx <- Sxtot_micromolperliter*Dliq*(SoilM/100)^(3.0) ## derived constant and parameters
  MMSx <- Sx/(KMSx+Sx)
  VmaxSx <- (alphaVmaxSx)*exp(-EaVmaxSx/(R*(SoilT+273.15))) 
  
  if(verbose>3) {
    cat("Dliq=",Dliq," ",
        "VmaxCH4prod=",VmaxCH4prod," ",
        "VmaxCH4ox=",VmaxCH4ox," ",
        "VmaxN2Oprod_nitrif=",VmaxN2Oprod_nitrif," ",
        "VmaxN2Oprod_denitrif=",VmaxN2Oprod_denitrif," ",
        "VmaxN2Ored=",VmaxN2Ored," ",
        "microsite.vol=",microsite.vol," ",
        "microsite.area=",microsite.area," ",
        "D_P=",D_P," ",
        "DO2=",DO2," ",
        "DCH4=",DCH4," ",
        "DO2=",DO2," ",
        "DN2O=",DN2O," ",
        "diffused_NO3=",diffused_NO3," ",
        "diffused_NH4=",diffused_NH4," ",
        "kl_CH4prod=",kl_CH4prod," ",
        "kMO2_Rh=",kMO2_Rh," ",
        "KMCH4ox=",KMCH4ox," ", ### changed from KMCH4ox_n to KMCH4ox
        "kMO2_CH4ox=",kMO2_CH4ox," ",
        "Sx=",Sx," ",
        "MMSx=",MMSx," ",
        "VmaxSx=",VmaxSx,"\n") ###Include Depth_O2 and Depth_CH4
  }
  
  ## Initialize

  #Need to change here
  iO2 <- O2airfrac*DO2*D_P #(for initialization)
  iO2_l_micromolperliter <- iO2*Henry_constant_O2*(10^6)
  iO2_l_micromolpermicrosite <- iO2_l_micromolperliter*(SoilM/100)*(1/1000)*microsite.vol
  iRh_l_micromolpermicrosite <- VmaxSx*MMSx*(iO2_l_micromolperliter/(KmO2_Rh+iO2_l_micromolperliter))*(SoilM/100)*(1/1000)*(microsite.vol)*(time.step)
  iRh <- iRh_l_micromolpermicrosite*(100/SoilM)*(10^3)*(1/microsite.vol)
  iO2_l_micromolpermicrosite_conc <- iO2_l_micromolpermicrosite-(iRh_l_micromolpermicrosite*(1/RQ))
  iO2_l_micromolperliter_conc <- iO2_l_micromolpermicrosite_conc*1000*(100/SoilM)*(1/microsite.vol)

  #Need to change here
  iCH4 <- CH4airfrac*DCH4*D_P #(for initialization)
  iCH4_l_micromolperliter <- iCH4*Henry_constant_CH4
  iCH4_l_micromolpermicrosite <- iCH4_l_micromolperliter*(SoilM/100)*(1/1000)*microsite.vol
  iCH4prod <- ((VmaxCH4prod*(1/1000)*(SoilM/100)*microsite.vol)/(1+(iO2_l_micromolperliter_conc/kl_CH4prod)))*(Sx/(Sx+KMSx_CH4prod))*(time.step)
  iCH4_l_conc_micromolpermicrosite <- iCH4_l_micromolpermicrosite+iCH4prod 
  iCH4_l_conc_micromolperliter <- iCH4_l_conc_micromolpermicrosite*1000*(100/SoilM)*(1/microsite.vol)
  iCH4ox <- ((VmaxCH4ox*(1/1000)*(SoilM/100)*microsite.vol))*(iCH4_l_conc_micromolperliter/(iCH4_l_conc_micromolperliter+KMCH4ox))*(iO2_l_micromolperliter_conc/(iO2_l_micromolperliter_conc+KmO2_CH4ox))*time.step 
  iCH4Flux_l_micromolpermicrosite <- iCH4prod-iCH4ox
  iCH4Flux_l_micromolperliter <- iCH4Flux_l_micromolpermicrosite*(100/SoilM)*(10^3)*(1/microsite.vol) 
  
  iN2O <- N2Oairfrac*DN2O*D_P #(for initialization)
  iN2O_l_micromolperliter <- iN2O*Henry_constant_N2O
  iN2O_l_micromolpermicrosite <- iN2O_l_micromolperliter*(SoilM/100)*(1/1000)*microsite.vol
  iN2Oprod_nitrif <- 0.5*((VmaxN2Oprod_nitrif*(1/1000)*(SoilM/100)*microsite.vol)/(1+(iO2_l_micromolperliter_conc/kl_N2Oprod)))*(diffused_NH4/(diffused_NH4+KM_N2Oprod_nitrif))*(iO2_l_micromolperliter_conc/(iO2_l_micromolperliter_conc+KMO2_N2Oprod_nitrif))*(time.step)
  iN2Oprod_denitrif <- 0.5*((VmaxN2Oprod_denitrif*(1/1000)*(SoilM/100)*microsite.vol)/(1+(iO2_l_micromolperliter_conc/kl_N2Oprod)))*(diffused_NO3/(diffused_NO3+KM_N2Oprod_denitrif))*(Sx/(Sx+KMSx_denitrif))*(time.step)
  iN2O_l_conc_micromolpermicrosite <- iN2O_l_micromolpermicrosite+iN2Oprod_nitrif+iN2Oprod_denitrif 
  iN2O_l_conc_micromolperliter <- iN2O_l_conc_micromolpermicrosite*1000*(100/SoilM)*(1/microsite.vol)
  iN2Ored <- ((VmaxN2Ored*(1/1000)*(SoilM/100)*microsite.vol)/(1+(iO2_l_micromolperliter_conc/kl_N2Ored)))*(iN2O_l_conc_micromolperliter/(iN2O_l_conc_micromolperliter+KMN2Ored))*(Sx/(Sx+KMSx_denitrif))*time.step 
  iN2OFlux_l_micromolpermicrosite <- iN2Oprod_nitrif+iN2Oprod_denitrif-iN2Ored
  iN2OFlux_l_micromolperliter <- iN2OFlux_l_micromolpermicrosite*(100/SoilM)*(10^3)*(1/microsite.vol) 
  
  deltaO2 <- 9999.0+reltol 
  deltaCH4 <- 9999.0+reltol 
  deltaN2O <- 9999.0+reltol 
  iter <- 0 
  
  #Need to change here
  if(verbose>0) cat("iter ",iter," O2=",iO2," O2_l_micromolperliter=",iO2_l_micromolperliter,
                    " O2_l_micromolpermicrosite=",iO2_l_micromolpermicrosite," 
                    Rh_l_micromolpermicrosite=",iRh_l_micromolpermicrosite," Rh=",iRh,
                    " O2_l_micromolpermicrosite_conc=",iO2_l_micromolpermicrosite_conc,
                    " O2_l_micromolperliter_conc=",iO2_l_micromolperliter_conc,
                    " CH4=",iCH4,"CH4_l_micromolperliter=",iCH4_l_micromolperliter,
                    "CH4_l_micromolpermicrosite=",iCH4_l_micromolpermicrosite," CH4prod=", iCH4prod,
                    " CH4_l_conc_micromolpermicrosite=", iCH4_l_conc_micromolpermicrosite,
                    " CH4_l_conc_micromolperliter=", iCH4_l_conc_micromolperliter,
                    " CH4ox=",iCH4ox,"CH4Flux_l_micromolpermicrosite=",iCH4Flux_l_micromolpermicrosite, 
                    "CH4Flux_l_micromolperliter=",iCH4Flux_l_micromolperliter,
                    
                    " N2O=",iN2O,"N2O_l_micromolperliter=",iN2O_l_micromolperliter,
                    "N2O_l_micromolpermicrosite=",iN2O_l_micromolpermicrosite,
                    "N2Oprod_nitrif=",iN2Oprod_nitrif,"N2Oprod_denitrif=",iN2Oprod_denitrif,
                    "N2O_l_conc_micromolpermicrosite=",iN2O_l_conc_micromolpermicrosite,
                    "N2O_l_conc_micromolperliter=",iN2O_l_conc_micromolperliter,
                    "N2Ored=",iN2Ored,"N2OFlux_l_micromolpermicrosite=",iN2OFlux_l_micromolpermicrosite, 
                    "N2OFlux_l_micromolperliter=",iN2OFlux_l_micromolperliter,"\n") 
  
  max_iter__ <- 99999
  max_iter_flag <- 0
  
  if(flagO2){
    ## If update O2, run O2 iteration until steady state
    ## run simulation until steady state 
    #while( (deltaCH4 > reltol) || (deltaO2 > reltol) )
    while( (deltaO2 > reltol))
    {
      if(iter > max_iter__){
        max_iter_flag <- 1
        break
      }

      #Need to change here
      O2_n_tmp <- O2airfrac*DO2*D_P #this is just to it identical to 1st iteration value
      O2 <- O2_n_tmp
      O2_l_micromolperliter <- O2*Henry_constant_O2*(10^6)
      O2_l_micromolpermicrosite <- O2_l_micromolperliter*(SoilM/100)*(1/1000)*microsite.vol
      Rh_l_micromolpermicrosite <- VmaxSx*MMSx*(O2_l_micromolperliter/(KmO2_Rh+O2_l_micromolperliter))*(SoilM/100)*(1/1000)*(microsite.vol)*(time.step)
      Rh <- Rh_l_micromolpermicrosite*(100/SoilM)*(10^3)*(1/microsite.vol)
      O2_l_micromolpermicrosite_conc <- O2_l_micromolpermicrosite-(Rh_l_micromolpermicrosite*(1/RQ))
      O2_l_micromolperliter_conc <- O2_l_micromolpermicrosite_conc*1000*(100/SoilM)*(1/microsite.vol)     
      
       ## Update delta
      deltaO2 <- abs(iO2-O2)/iO2 
      ## update states

      iO2 <- O2
      iO2_l_micromolperliter <- O2_l_micromolperliter
      iO2_l_micromolpermicrosite <- O2_l_micromolpermicrosite
      iRh_l_micromolpermicrosite <- Rh_l_micromolpermicrosite
      iRh <- Rh
      iO2_l_micromolpermicrosite_conc <- O2_l_micromolpermicrosite_conc 
      iO2_l_micromolperliter_conc <- O2_l_micromolperliter_conc

      iter <- iter + 1 
      if(verbose>2) cat("iter ",iter," O2=",iO2," O2_l_micromolperliter=",iO2_l_micromolperliter,
                        " O2_l_micromolpermicrosite=",iO2_l_micromolpermicrosite,   
                        " Rh_l_micromolpermicrosite=",iRh_l_micromolpermicrosite," Rh=",iRh,
                        " O2_l_micromolpermicrosite_conc=",iO2_l_micromolpermicrosite_conc,
                        " O2_l_micromolperliter_conc=",iO2_l_micromolperliter_conc,
                        "(delta=",deltaO2,")","\n") 
    } ## End O2 iterations
  
  } ## End update O2
  else{
    ## keep initial values
    O2 <- iO2
    O2_l_micromolperliter <- iO2_l_micromolperliter
    O2_l_micromolpermicrosite <-  iO2_l_micromolpermicrosite
    Rh_l_micromolpermicrosite <- iRh_l_micromolpermicrosite
    Rh <- iRh
    O2_l_micromolpermicrosite_conc <- iO2_l_micromolpermicrosite_conc 
    O2_l_micromolperliter_conc <- iO2_l_micromolperliter_conc

  } ## End not update O2
  if(verbose>0) cat("###### End O2 Iterations ##### \n iter = ",
                    iter," O2=",iO2,"(delta=",deltaO2,")","\n") 
  iter <- 0  ## re-set iterations for CH4
  
  if(flagCH4){
    ## If update CH4, run CH4 iteration until steady state
    ## run simulation until steady state 
    #while( (deltaCH4 > reltol) || (deltaO2 > reltol) )
    while( (deltaCH4 > reltol))
    {
      if(iter > max_iter__){
        max_iter_flag <- 1
        break
      }

      #Need to change here
      CH4 <- CH4airfrac*DCH4*D_P #this is just to it identical to 1st iteration value
      CH4_l_micromolperliter <- CH4*Henry_constant_CH4
      CH4_l_micromolpermicrosite <- CH4_l_micromolperliter*(SoilM/100)*(1/1000)*microsite.vol
      CH4prod <- ((VmaxCH4prod*(1/1000)*(SoilM/100)*microsite.vol)/(1+(O2_l_micromolperliter_conc/kl_CH4prod)))*(Sx/(Sx+KMSx_CH4prod))*(time.step)
      CH4_l_conc_micromolpermicrosite <- CH4_l_micromolpermicrosite+CH4prod 
      CH4_l_conc_micromolperliter <- CH4_l_conc_micromolpermicrosite*1000*(100/SoilM)*(1/microsite.vol)
      CH4ox <- ((VmaxCH4ox*(1/1000)*(SoilM/100)*microsite.vol))*(CH4_l_conc_micromolperliter/(CH4_l_conc_micromolperliter+KMCH4ox))*(O2_l_micromolperliter_conc/(O2_l_micromolperliter_conc+KmO2_CH4ox))*time.step 
      CH4Flux_l_micromolpermicrosite <- CH4prod-CH4ox
      CH4Flux_l_micromolperliter <- CH4Flux_l_micromolpermicrosite*(100/SoilM)*(10^3)*(1/microsite.vol) 
      

      ##### Need to change here to end
      
      ## Update delta
      deltaCH4 <- abs(iCH4-CH4)/iCH4 

      iCH4 <- CH4
      iCH4_l_micromolperliter <- CH4_l_micromolperliter
      iCH4_l_micromolpermicrosite <- CH4_l_micromolpermicrosite 
      iCH4prod <- CH4prod
      iCH4_l_conc_micromolpermicrosite <- CH4_l_conc_micromolpermicrosite
      iCH4_l_conc_micromolperliter <- CH4_l_conc_micromolperliter
      iCH4ox <- CH4ox
      iCH4Flux_l_micromolpermicrosite <- CH4Flux_l_micromolpermicrosite
      iCH4Flux_l_micromolperliter <- CH4Flux_l_micromolperliter
        
      iter <- iter + 1 
      if(verbose>2) cat("iter ",iter," CH4=",iCH4," CH4_l_micromolperliter=", iCH4_l_micromolperliter,
                        " CH4_l_micromolpermicrosite=", iCH4_l_micromolpermicrosite," CH4prod=", iCH4prod,
                        " CH4_l_conc_micromolpermicrosite=", iCH4_l_conc_micromolpermicrosite,
                        " CH4_l_conc_micromolperliter=", iCH4_l_conc_micromolperliter,
                        " CH4ox=",iCH4ox,"CH4Flux_l_micromolpermicrosite=",iCH4Flux_l_micromolpermicrosite,
                        "CH4Flux_l_micromolperliter=",iCH4Flux_l_micromolperliter,
                        "(deltaCH4=",deltaCH4,")\n") 
    } ## End iterations
   
  } ## End updating CH4
  else{
    ## keep initial values
    CH4 <- iCH4
    CH4_l_micromolperliter <- iCH4_l_micromolperliter
    CH4_l_micromolpermicrosite <- iCH4_l_micromolpermicrosite 
    CH4prod <- iCH4prod
    CH4_l_conc_micromolpermicrosite <- iCH4_l_conc_micromolpermicrosite
    CH4_l_conc_micromolperliter <- iCH4_l_conc_micromolperliter
    CH4ox <- iCH4ox
    CH4Flux_l_micromolpermicrosite <- iCH4Flux_l_micromolpermicrosite
    CH4Flux_l_micromolperliter <- iCH4Flux_l_micromolperliter    
    
  } ## End not updating CH4
  if(verbose>0) cat("###### End CH4 Iterations ##### \n iter = ",
                    iter," CH4=",iCH4,"(delta=",deltaCH4,")","\n") 
  iter <- 0  ## re-set iterations for N20
  
  if(flagN2O){
    ## If update N2O, run N2O iteration until steady state
    ## run simulation until steady state
    #while( (deltaN2O > reltol) || (deltaO2 > reltol) )
    while( (deltaN2O > reltol))
    {
      if(iter > max_iter__){
        max_iter_flag <- 1
        break
      }

      N2O <- N2Oairfrac*DN2O*D_P #this is identical to 1st iteration
      N2O_l_micromolperliter <- N2O*Henry_constant_N2O
      N2O_l_micromolpermicrosite <- N2O_l_micromolperliter*(SoilM/100)*(1/1000)*microsite.vol
      N2Oprod_nitrif <- 0.5*((VmaxN2Oprod_nitrif*(1/1000)*(SoilM/100)*microsite.vol)/(1+(O2_l_micromolperliter_conc/kl_N2Oprod)))*(diffused_NH4/(diffused_NH4+KM_N2Oprod_nitrif))*(O2_l_micromolperliter_conc/(O2_l_micromolperliter_conc+KMO2_N2Oprod_nitrif))*(time.step)
      N2Oprod_denitrif <- 0.5*((VmaxN2Oprod_denitrif*(1/1000)*(SoilM/100)*microsite.vol)/(1+(O2_l_micromolperliter_conc/kl_N2Oprod)))*(diffused_NO3/(diffused_NO3+KM_N2Oprod_denitrif))*(Sx/(Sx+KMSx_denitrif))*(time.step)
      N2O_l_conc_micromolpermicrosite <- N2O_l_micromolpermicrosite+N2Oprod_nitrif+N2Oprod_denitrif
      N2O_l_conc_micromolperliter <- N2O_l_conc_micromolpermicrosite*1000*(100/SoilM)*(1/microsite.vol)
      N2Ored <- ((VmaxN2Ored*(1/1000)*(SoilM/100)*microsite.vol)/(1+(O2_l_micromolperliter_conc/kl_N2Ored)))*(N2O_l_conc_micromolperliter/(N2O_l_conc_micromolperliter+KMN2Ored))*(Sx/(Sx+KMSx_denitrif))*time.step
      N2OFlux_l_micromolpermicrosite <- N2Oprod_nitrif+N2Oprod_denitrif-N2Ored
      N2OFlux_l_micromolperliter <- N2OFlux_l_micromolpermicrosite*(100/SoilM)*(10^3)*(1/microsite.vol)

      ## Update delta
      deltaN2O <- abs(iN2O-N2O)/iN2O

      iN2O <- N2O
      iN2O_l_micromolperliter <- N2O_l_micromolperliter
      iN2O_l_micromolpermicrosite <- N2O_l_micromolpermicrosite
      iN2Oprod_nitrif <- N2Oprod_nitrif
      iN2Oprod_denitrif <- N2Oprod_denitrif
      iN2O_l_conc_micromolpermicrosite <- N2O_l_conc_micromolpermicrosite
      iN2O_l_conc_micromolperliter <- N2O_l_conc_micromolperliter
      iN2Ored <- N2Ored
      iN2OFlux_l_micromolpermicrosite <- N2OFlux_l_micromolpermicrosite
      iN2OFlux_l_micromolperliter <- N2OFlux_l_micromolperliter

      iter <- iter + 1
      if(verbose>2) cat("iter ",iter," N2O=",iN2O," N2O_l_micromolperliter=", iN2O_l_micromolperliter,
                        " N2O_l_micromolpermicrosite=", iN2O_l_micromolpermicrosite,
                        "N2Oprod_nitrif=",iN2Oprod_nitrif,"N2Oprod_denitrif=",iN2Oprod_denitrif,
                        " N2O_l_conc_micromolpermicrosite=", iN2O_l_conc_micromolpermicrosite,
                        " N2O_l_conc_micromolperliter=", iN2O_l_conc_micromolperliter,
                        " N2Ored=",iN2Ored,"N2OFlux_l_micromolpermicrosite=",iN2OFlux_l_micromolpermicrosite,
                        "N2OFlux_l_micromolperliter=",iN2OFlux_l_micromolperliter,
                        "(deltaN2O=",deltaN2O,")\n")
    } ## End iterations
  } ## End updating N2O
  else{
    ## keep initial values
    N2O <- iN2O
    N2O_l_micromolperliter <- iN2O_l_micromolperliter
    N2O_l_micromolpermicrosite <- iN2O_l_micromolpermicrosite
    N2Oprod_nitrif <- iN2Oprod_nitrif
    N2Oprod_denitrif <- iN2Oprod_denitrif
    N2O_l_conc_micromolpermicrosite <- iN2O_l_conc_micromolpermicrosite
    N2O_l_conc_micromolperliter <- iN2O_l_conc_micromolperliter
    N2Ored <- iN2Ored
    N2OFlux_l_micromolpermicrosite <- iN2OFlux_l_micromolpermicrosite
    N2OFlux_l_micromolperliter <- iN2OFlux_l_micromolperliter

  } ## End not updating N2O
  if(verbose>0) cat("###### End N2O Iterations ##### \n iter = ",
                    iter," N2O=",iN2O,"(delta=",deltaN2O,")","\n")
  ## Return steady state
  damm_return <- rep(NA,28)
 
  
  names(damm_return) <- c("O2","O2_l_micromolperliter","O2_l_micromolpermicrosite","Rh_l_micromolpermicrosite","Rh",
                        "O2_l_micromolpermicrosite_conc","O2_l_micromolperliter_conc","CH4","CH4_l_micromolperliter",
                      "CH4_l_micromolpermicrosite","CH4prod","CH4_l_conc_micromolpermicrosite",
                     "CH4_l_conc_micromolperliter","CH4ox","CH4Flux_l_micromolpermicrosite",
                      "CH4Flux_l_micromolperliter","D_P","flag",
                     "N2O","N2O_l_micromolperliter",
                     "N2O_l_micromolpermicrosite","N2Oprod_nitrif","N2Oprod_denitrif","N2O_l_conc_micromolpermicrosite",
                     "N2O_l_conc_micromolperliter","N2Ored","N2OFlux_l_micromolpermicrosite",
                     "N2OFlux_l_micromolperliter")
  damm_return[1] <- O2 
  damm_return[2] <- O2_l_micromolperliter  
  damm_return[3] <- O2_l_micromolpermicrosite 
  damm_return[4] <- Rh_l_micromolpermicrosite
  damm_return[5] <- Rh
  damm_return[6] <- O2_l_micromolpermicrosite_conc
  damm_return[7] <- O2_l_micromolperliter_conc
  damm_return[8] <- CH4 
  damm_return[9] <- CH4_l_micromolperliter 
  damm_return[10] <- CH4_l_micromolpermicrosite 
  damm_return[11] <- CH4prod 
  damm_return[12] <- CH4_l_conc_micromolpermicrosite
  damm_return[13] <- CH4_l_conc_micromolperliter
  damm_return[14] <- CH4ox
  damm_return[15] <- CH4Flux_l_micromolpermicrosite
  damm_return[16] <- CH4Flux_l_micromolperliter
  damm_return[17] <- D_P
  damm_return[18] <- max_iter_flag
  
  damm_return[19] <- N2O 
  damm_return[20] <- N2O_l_micromolperliter 
  damm_return[21] <- N2O_l_micromolpermicrosite 
  damm_return[22] <- N2Oprod_nitrif 
  damm_return[23] <- N2Oprod_denitrif
  damm_return[24] <- N2O_l_conc_micromolpermicrosite
  damm_return[25] <- N2O_l_conc_micromolperliter
  damm_return[26] <- N2Ored
  damm_return[27] <- N2OFlux_l_micromolpermicrosite
  damm_return[28] <- N2OFlux_l_micromolperliter
  return(damm_return) 
}
test1 <- function()
{
  rm(list=ls()) 
  ## Session > Work Dir > Source File Loc
  source("damm_11a-non-steady-state.R")
  microsite <- read.csv("microsite_5.csv")
  sink("./damm_liq_debug.txt")
  r <- damm_work( 
    R=0.0083, O2airfrac=0.209,  CH4airfrac=1.8,N2Oairfrac=0.328,BD=1.009, PD=2.52,  Porosity=0.5220,
    p=0.00046, RQ=1, 
    
    SoilM.o=30,   SoilT=20,   Sxtot=0.0010, SoilM=22.5, 
    Depth_O2=10,Depth_CH4=10,Depth_N2O=10,#Depth ,  ###correct for diff depth
    
    total.microsite=10000,NO3=0.7, NH4=5,
    time.step=300,
    
    iO2=0.209, iCH4=1.8, iN2O=0.328, ## initial O2 and CH4 values
    
    flagO2=TRUE,  flagCH4=TRUE, flagN2O=TRUE,## whether run O2 and CH4 to steady states
    
    KmO2_Rh=245, alphaVmaxSx=230000000000,   EaVmaxSx=72,   KMSx=1000,   
    alphaVmaxCH4prod=3000000, EaVmaxCH4prod=100,  KMSx_CH4prod=350,
    alphaVmaxCH4ox=0.07, EaVmaxCH4ox=30,   KMCH4ox=0.01, KmO2_CH4ox=43, kl_CH4prod=3,
    
    alphaVmaxN2Oprod_nitrif=100, EaVmaxN2Oprod_nitrif=50,  KM_N2Oprod_nitrif=22.22,KMO2_N2Oprod_nitrif=143.3,
    alphaVmaxN2Oprod_denitrif=100, EaVmaxN2Oprod_denitrif=50,  KM_N2Oprod_denitrif=26,
    alphaVmaxN2Ored=1500, EaVmaxN2Ored=50,   KMN2Ored=0.16, 
    kl_N2Oprod=14.3,kl_N2Ored=7.2,KMSx_denitrif=1000,
    
    reltol=sqrt(.Machine$double.eps),verbose=99
  )
  sink()
  r
 }

#for test
R=0.0083 
O2airfrac=0.209 
CH4airfrac=1.8
N2Oairfrac=0.328
BD=1.009
PD=2.52
Porosity=0.5220
p=0.00046
RQ=1
SoilM.o=30
SoilT=20
Sxtot=0.0010
SoilM=22.5
Depth_O2=10
Depth_CH4=10
Depth_N2O=10#Depth ,  ###correct for diff depth

total.microsite=10000
NO3=0.7
NH4=5
time.step=300

iO2=0.209
iCH4=1.8
iN2O=0.328 ## initial O2 and CH4 values

flagO2=TRUE; flagCH4=TRUE;flagN2O=TRUE## whether run O2 and CH4 to steady states

KmO2_Rh=245; alphaVmaxSx=230000000000;   EaVmaxSx=72;   KMSx=1000   
alphaVmaxCH4prod=3000000; EaVmaxCH4prod=100;  KMSx_CH4prod=350
alphaVmaxCH4ox=0.07;EaVmaxCH4ox=30;   KMCH4ox=0.01; KmO2_CH4ox=43; kl_CH4prod=3

alphaVmaxN2Oprod_nitrif=100; EaVmaxN2Oprod_nitrif=50;  KM_N2Oprod_nitrif=22.22;KMO2_N2Oprod_nitrif=143.3
alphaVmaxN2Oprod_denitrif=100; EaVmaxN2Oprod_denitrif=50;  KM_N2Oprod_denitrif=26
alphaVmaxN2Ored=1500; EaVmaxN2Ored=50;   KMN2Ored=0.16 
kl_N2Oprod=14.3;kl_N2Ored=7.2;KMSx_denitrif=1000

reltol=sqrt(.Machine$double.eps);verbose=99