#This script runs the carbon flux and isotope model based on 'yourdata', which could be empirical
#culture data, or generic artificial data to investigate model behavior.

######################################
##Required packages
######################################
packages <- c("seacarb","readxl","ggplot2","dplyr")
lapply(packages, require, character.only = TRUE)

######################################
##grind.R, to solve ODEs
######################################
source("grind.R") #download and save to working directory: https://bioinformatics.bio.uu.nl/rdb/grind.html

######################
##Model loop
######################
solver_model <- function(yourdata){

results_list <- vector("list", nrow(yourdata))  

for(i in c(1:nrow(yourdata))){
  print(paste("Processing row:", i))
  
  #######################
  ##Parameters
  #######################
  source("parameters.R",local=TRUE) #External file for parameter settings
  
  ######################
  ##Carbon flux model
  ######################
  carbon_model <- function(t, state, parms) {
    with(as.list(c(state,parms)), { 
  
      
      dCc <- kf1 * Hc - kr1 * Cc + SAcyt/Vcyt*aE*(Ce)-SAcyt/Vcyt*aE*(Cc) - SAch/Vcyt*aC*(Cc)+SAch/Vcyt*aC*(Cch) + 0.5*0.7*calc*Hc
      dCch <-kf2 * Hch - kr2 * Cch + SAch/Vch*aC * (Cc)-SAch/Vch*aC * (Cch) - SAp/Vch*pP*(Cch)+SAp/Vch*pP*(Cp)
      dHc <- kr1*Cc - kf1*Hc - uC1*(Hc)/(Kn+(Hc))*SAch/Vcyt + uE*He/(Kn+He)*SAcyt/Vcyt - 0.85*calc*Hc
      dHch <- kr2*Cch - kf2*Hch + uC1*(Hc)/(Kn+(Hc))*SAch/Vch - uC2*(Hch)/(Kn+(Hch))*SAt/Vch
      
      dCt<-kf3*Ht-kr3*Ct - SAp/Vt*pT*(Ct)+SAp/Vt*pT*(Cp)
      dHt<-kr3*Ct-kf3*Ht+uC2*(Hch)/(Kn+(Hch))*SAt/Vt
      dCp<-SAp/Vp*pP*(Cch)-SAp/Vp*pP*(Cp) + SAp/Vp*pT*(Ct)-SAp/Vp*pT*(Cp) - r1 + kf4*Hp-kr4*Cp
      dHp<-kr4*Cp-kf4*Hp
      
      dCaCO3<-0.85*calc*Hc
      
      return(list(c(dCc,dCch,dHc,dHch,dCt,dHt,dCp,dHp,dCaCO3)))
    }) 
  }  
  
  p <- c(calc=calc,
         r1=r1,
         aE=pC,aC=pC,pT=pC,pP=pP,r_ch=r_ch,
         r_cyt=r_cyt,Kn=Kn,uC1=uptake.rate,uC2=uptake.rate,
         uE=uptake.rate,Ce=Ce,He=He,Temp=Temp)
  s<-c(Cc=Ce,Cch=Ce,Hc=He,Hch=He,Ct=0,Ht=0,Cp=0,Hp=0,CaCO3=0)
  
  r<-run(odes=carbon_model,state = s,parms = p,1e8,tstep=1e7,timeplot = F)
  
  ############################
  ## Steady state conditions
  ############################
  Hc <- as.numeric(r["Hc"]) 
  Cc <- as.numeric(r["Cc"])     
  Cch <-as.numeric(r["Cch"])   
  Hch <- as.numeric(r["Hch"])    
  Ct<-as.numeric(r["Ct"])
  Ht<-as.numeric(r["Ht"])
  Cp<-as.numeric(r["Cp"])
  Hp<-as.numeric(r["Hp"])
  
  #########################
  ##Determining relative HCO3- uptake and CO2 leakage
  #########################
  uE=as.numeric(p["uE"])
  diff<-(SAcyt/Vcyt*pC*(Ce))
  eff<-SAcyt/Vcyt*pC*(Cc)
  upt<- uE*He/(Kn+He)*SAcyt/Vcyt
  Uptake<-(upt/((diff)+upt))
  leakage<-eff/(diff+eff)
  
  ########################
  ##Isotope model
  ########################
  isotope_model <- function(t, state, parms) {
    with(as.list(c(state,parms)), {
      
      ddCc <- kf1 * Hc * (dHc -fracHC) - kr1 * Cc * (dCc - fracCH) + SAcyt/Vcyt*aE*(Ce)*dCe-SAcyt/Vcyt*aE*(Cc)*dCc - SAch/Vcyt*aC*(Cc)*dCc+SAch/Vcyt*aC*(Cch)*dCch + 0.5*0.7*calc*Hc*(dHc-1)
      ddCch <-kf2 * Hch*(dHch-fracHC) - kr2 * Cch*(dCch-fracCH) + SAch/Vch*aC * (Cc)*dCc-SAch/Vch*aC * (Cch)*dCch - SAp/Vch*pP*(Cch)*dCch+SAp/Vch*pP*(Cp)*dCp 
      ddHc <- kr1*Cc*(dCc - fracCH) - kf1*Hc*(dHc -fracHC) - uC1*dHc*Hc/(Kn+Hc)*SAch/Vcyt + uE*dHe*He/(Kn+He)*SAcyt/Vcyt - 0.85*calc*Hc*(dHc+1)
      ddHch <- kr2*Cch*(dCch-fracCH) - kf2*Hch*(dHch-fracHC) + uC1*dHc*Hc/(Kn+Hc)*SAch/Vch - uC2*Hch/(Kn+Hch)*SAt/Vch*dHch
      
      ddCt<-kf3*Ht*(dHt-fracHC)-kr3*Ct*(dCt-fracCH) - SAp/Vt*pT*(Ct)*dCt+SAp/Vt*pT*(Cp)*dCp 
      ddHt<-kr3*Ct*(dCt-fracCH)-kf3*Ht*(dHt-fracHC)+uC2*Hch/(Kn+Hch)*SAt/Vt*dHch
      ddCp<-SAp/Vp*pP*(Cch)*dCch-SAp/Vp*pP*(Cp)*dCp + SAp/Vp*pT*(Ct)*dCt-SAp/Vp*pT*(Cp)*dCp- r1*(dCp-rubp.frac)+kf4*Hp*(dHp-fracHC)-kr4*Cp*(dCp-fracCH)
      ddHp<-kr4*Cp*(dCp-fracCH) - kf4*Hp*(dHp-fracHC)
      
      ddCaCO3<-0.85*calc*Hc*(dHc-dCaCO3+1)
      
      return(list(c(ddCc,ddCch,ddHc,ddHch,ddCt,ddHt,ddCp,ddHp,ddCaCO3)))
    })
  }
  p <- c(calc=calc,Hc = as.numeric(r["Hc"]),Cc = as.numeric(r["Cc"]),Cch = as.numeric(r["Cch"])  ,Hch = as.numeric(r["Hch"]),
         kf1=kf1,kr1=kr1,kf2=kf2,kr2=kr2,pC=pC,r1=r1,kf3=kf3,kr3=kr3,kf4=kf4,kr4=kr4,
         aE=pC,aC=pC,r_ch=r_ch,r_cyt=r_cyt,Kn=Kn,uC1=uptake.rate,uC2=uptake.rate,pT=pC,
         uC3=uptake.rate,uE=uptake.rate,Ce=Ce,He=He,dHe=dHe,dCe=dCe)
  s<-c(dCc=0,dCch=0,dHc=0,dHch=0,dCt=0,dHt=0,dCp=0,dHp=0,dCaCO3=0)
  r2<-run(odes=isotope_model,state=s,parms=p,tmax=2e14,tstep=1e13,ymax=40,ymin=-20,timeplot = F)
  
  ##############################
  ##Compartment-specific d13C of CO2 and HCO3-
  ##############################
  dCc<-as.numeric(r2["dCc"])
  dCch<-as.numeric(r2["dCch"])
  dHc<-as.numeric(r2["dHc"])
  dHch<-as.numeric(r2["dHch"])
  dCp<-as.numeric(r2["dCp"])
  dHp<-as.numeric(r2["dHp"])
  dHt<-as.numeric(r2["dHt"])
  dCt<-as.numeric(r2["dCt"])
  dCaCO3<-as.numeric(r2["dCaCO3"])
  EpcalcCO2<-dCe-dCaCO3
  
  ##########################
  ##Ep with sugar offset
  #########################
  eps_eff<-rubp.frac+2

  dB<-dCp-eps_eff
  dBtot<-dB
  eps=(dCe-dBtot)
  
  ref=yourdata$ref[i]
  species=yourdata$species[i]
  
  #######################
  ##Save results in data frame
  #######################
  if('ref' %in% colnames(yourdata)){
  results_list[[i]] <- data.frame(
    dCc=dCc, dCch=dCch, dHc=dHc, dHch=dHch, dCt=dCt, dCp=dCp, dHt=dHt, dHp=dHp, 
    He=He, dBtot=dBtot, Ce=Ce, pCO2=pCO2, Uptake=Uptake, Temp=Temp, 
    daylight=daylight, eps=eps, PFD=PFD, r1=r1, Leakage=leakage, 
    Vcell=Vcell, Cc=Cc, Cch=Cch, Ct=Ct, Cp=Cp, species=species, 
    dCaCO3=dCaCO3, EpcalcCO2=EpcalcCO2, ref=ref
  )}else{
    results_list[[i]] <- data.frame(
      dCc=dCc, dCch=dCch, dHc=dHc, dHch=dHch, dCt=dCt, dCp=dCp, dHt=dHt, dHp=dHp, 
      He=He, dBtot=dBtot, Ce=Ce, pCO2=pCO2, Uptake=Uptake, Temp=Temp, 
      daylight=daylight, eps=eps, PFD=PFD, r1=r1, Leakage=leakage, 
      Vcell=Vcell, Cc=Cc, Cch=Cch, Ct=Ct, Cp=Cp, species=species, 
      dCaCO3=dCaCO3, EpcalcCO2=EpcalcCO2)
  }
  
  
}
  df_model <- do.call(rbind, results_list)
  return(df_model)
}
