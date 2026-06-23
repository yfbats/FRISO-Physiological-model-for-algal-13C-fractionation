#This file contains all the parameters used to run the Carbon flux and isotope model.
#If using culture data (culture.data) to fit/compare model to, the constrained parameters are 
#extracted and determined, whereas the tunable parameters (line 114; calc, uptake.rate, pPfactor) can be set freely.

#If studying the generic behavior of the model, all parameters can be set freely.
#Generic artificial data for running model should be stored to 'yourdata'.

#######################
##Rubisco fractionation
#######################
if(yourdata$species[i]=="E. huxleyi"){rubp.frac=11.1}else{rubp.frac=24}

#######################
##Membrane permeability (/cm)
#######################
if(yourdata$species[i]=="E. huxleyi"){pC=0.01}else{pC=0.027}

#######################
##Cell geometry
#######################
Vcell=yourdata$Vcell[i]

#Cytosol radius
r_cyt=((3*Vcell*0.847)/(4*pi))^(1/3)

#Chloroplast radius
r_ch=0.47*r_cyt

#Pyrenoid radius
r_p=0.15*r_cyt

#Thylakoid radius
r_t=0.41*r_cyt

#Cytosol surface area
SAcyt=4*pi*r_cyt^2

#Chloroplast surface area
SAch=4*pi*r_ch^2

#Cytosol volume
Vcyt=4/3*pi*r_cyt^3

#Chloroplast volume
Vch=4/3*pi*r_ch^3

#Thylakoid surface area
SAt=4*pi*r_t^2

#Pyrenoid surface area
SAp=4*pi*r_p^2

#Thylakoid volume
Vt=4/3*pi*r_t^3

#Pyrenoid volume
Vp=4/3*pi*r_p^3


#######################
##Carbon fixation rate
#######################
if ("r1" %in% colnames(yourdata)) {
   r1=yourdata$r1[i]
  }else{
poc=yourdata$POC[i]
ui=as.numeric(yourdata$ui[i])
r1<-poc/Vcell*(ui/24/60/60)
  }

#######################
##Ambient (culture/environmental) conditions
#######################
PFD<-as.numeric(yourdata$PFD[i])
daylight=yourdata$daylight[i]
Temp=yourdata$Temp[i]
sal=yourdata$sal[i]
pH=yourdata$pH[i]

######################
##Compartment specific pH and reaction rates
######################
Ka=as.numeric(K1(T=Temp,S=sal)*1.025)#dissociation constant

pHcyt<- 7
pHch=8
pHt=6
pHp=8

kf1=10^(-pHcyt)/Ka*1e6
kr1=Ka/(10^(-pHcyt))*1e6
kf2=10^(-pHch)/Ka*1e6
kr2=Ka/10^(-pHch)*1e6
kf3=10^(-pHt)/Ka*1e6
kr3=Ka/10^(-pHt)*1e6
kf4=10^(-pHp)/Ka*1e6
kr4=Ka/10^(-pHp)*1e6

#######################
##Carbonate chemistry
#######################
Ce = yourdata$co2[i]/1e6*1.025/1000  # External CO2 concentration
He<-as.numeric(carb(flag=1,var2=Ce*1000/1.025,var1=pH,T=Temp,S=sal)[,14])*1.025/1000 #External HCO3- concentration
pCO2<-carb(flag=2,var1=Ce*1e3/1.025,var2=He*1e3/1.025,S=sal,T=Temp)[,8]

###Equilibrium fractionation between HCO3- and CO2 (aq)
eq.frac=9.866*(10^3/(Temp+273.15))-24.12 #(Mook et al 1974)
fracHC=eq.frac
fracCH=fracHC-eq.frac

###ambient d13C of HCO3- and CO2
dHe=0
dCe=dHe-eq.frac


#######################
##Tunable parameters
#######################
if(yourdata$species[i]=="E. huxleyi"){
  if(daylight==24){pPfactor=40}
  else{pPfactor=1}
  calc=0.07*daylight+0.007*PFD-1.3
  uptake.rate=1.1e-11
  Kn=2e-7

}else{
  
  if(daylight==24){pPfactor=10}
  else{pPfactor=1}
  calc=0
  Kn=0
  if(yourdata$species[i]=="Alexandrium"){uptake.rate=2e-10}
  if(yourdata$species[i]=="Gonyaulax"){uptake.rate=5e-10}
  if(yourdata$species[i]=="Scrippsiella"){uptake.rate=1e-9}
  if(yourdata$species[i]=="Protoceratium"){uptake.rate=9e-13*PFD+5e-10}
  
}
pP=pPfactor*pC

