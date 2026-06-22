#This script compares model output with culture data of G(E). huxleyi and the dinoflagellates
#A. tamarense, S. trochoidea, G. spinifera, and P. reticulatum.

###################################
##Set correct working directory
###################################
#setwd("")

######################################
##Required packages
######################################
packages <- c("seacarb","readxl","ggplot2","dplyr")
lapply(packages, require, character.only = TRUE)

######################################
##grind.R, to solve ODEs
######################################
source("grind.R") #download and save to working directory: https://bioinformatics.bio.uu.nl/rdb/grind.html

###########################
##Load in culture data
###########################

#Data should contain the following columns:
# daylight (daylength in h)
# co2 (in µmol/L)
# pH
# Vcell (volume of cell in cm^3)
# eps (ε[p] (fractionation) in per mil)
# Temp (temperature in degrees Celsius)
# sal (salinity)
# PFD (photon flux density, in µmol photons/m^2/s)
# POC (cellular carbon content in mol)
# ui (instantaneous growth rate in /d)
# species (Species names)
# EpcalcCO2 (in per mil)

#######################################################################
###Hash/unhash to select G(E). huxleyi or Dinoflagellate culture data:
#######################################################################

###G. huxleyi
culture.data<-as.data.frame(read_excel("data/Ehux.xlsx"))
fluxdata <- read_excel("data/haptfluxdata.xlsx")


###Dinoflagellates
culture.data<-as.data.frame(read_excel("data/Dino.xlsx"))
fluxdata <- read_excel("data/dinofluxdata.xlsx")

########################################
##run model and save results in df_model
########################################
source("Carbon flux and isotope model.R")
df_model=solver_model(culture.data)

##############################
##Check fit
##############################

df_comb<-data.frame(culture.data$eps,df_model$eps,culture.data$daylight,culture.data$co2,df_model$Ce*1e9/1.025,culture.data$ref)
colnames(df_comb)<-c("measured","modeled","daylight","co2meas","co2mod",'ref')

###Root mean square error###
rmse=sqrt(mean((df_comb$measured - df_comb$modeled)^2));rmse

###R squared###
resid2<-with(df_comb,(1-sum((measured-modeled)^2)/(sum((measured-mean(measured))^2))));resid2

ggplot(df_comb)+
  theme_bw()+
  geom_point(aes(x=measured,y=modeled),size=2)+
  geom_abline(slope=1,intercept=0,linetype='dashed')+
  xlab(expression(Measured~epsilon[p]~"(‰)"))+
  ylab(expression(Modeled~epsilon[p]~"(‰)"))+
  xlim(0,30)+
  ylim(0,30)+
  geom_text(aes(x=10,y=23),label=paste("RMSE = ",round(rmse,2),"‰"),size=7)+
  theme(axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        panel.grid = element_blank(),
        text = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.position = c(0.85, 0.25),  # inside plot (x,y)
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA)
  )

###################################
##Plot model versus culture data
###################################

#Bin PFD
df_model <- df_model %>%
  mutate(
    PFD = as.numeric(as.character(PFD)),   # convert factor → numeric
    PFD_group = case_when(
      PFD == 15 ~ "15",
      PFD >= 30 & PFD <= 45 ~ "30-45",
      PFD >= 80 & PFD <= 85 ~ "80-85",
      PFD >= 150 & PFD <= 160 ~ "150-160",
      PFD >= 175 & PFD <= 190 ~ "175-190",
      PFD == 250 ~ "250",
      TRUE ~ NA_character_
    ),
    PFD_group = factor(PFD_group, levels = rev(c("15", "30-45", "80-85", "150-160", "175-190","250")))
  )

culture.data <- culture.data %>%
  mutate(
    PFD = as.numeric(as.character(PFD)),
    PFD_group = case_when(
      PFD == 15 ~ "15",
      PFD >= 30 & PFD <= 45 ~ "30-45",
      PFD >= 80 & PFD <= 85 ~ "80-85",
      PFD >= 150 & PFD <= 160 ~ "150-160",
      PFD >= 175 & PFD <= 190 ~ "175-190",
      PFD == 250 ~ "250",
      TRUE ~ NA_character_
    ),
    PFD_group = factor(PFD_group, levels = rev(c("15", "30-45", "80-85", "150-160", "175-190","250")))
  )

if(unique(culture.data$species=="E. huxleyi")){
ggplot()+
  theme_bw()+
  ylim(0,26.7)+
  xlim(0,90)+
  geom_point(data=culture.data,aes(x=(co2*1.025),y=eps,color=(PFD_group),shape='Culture'),size=3)+
  geom_point(data=df_model[df_model$daylight==16,],aes(x=Ce*1e9,y=eps,color=(PFD_group),shape='Model'),size=4)+
  geom_point(data=df_model[df_model$daylight==24,],aes(x=Ce*1e9,y=eps,color=(PFD_group),shape='Model'),size=4)+
  geom_smooth(data=df_model[df_model$daylight==24,],method='lm',aes(x=Ce*1e9,y=eps,color=PFD_group),formula=y~poly(x,2),se=F,show.legend = F)+
  geom_smooth(data=df_model[df_model$daylight==16,],method='lm',aes(x=Ce*1e9,y=eps,color=PFD_group),formula=y~poly(x,2),se=F,show.legend = F)+
  xlab(expression(CO[2~"("*aq*")"]*" (μmol/L)"))+
  ylab(expression(epsilon[p]~"(‰)"))+
  
  geom_text(aes(x=7.2,y=17),label=expression("{"),size=18)+
  geom_text(aes(x=0.4,y=17),label=expression("24 h"),size=7)+
  geom_text(aes(x=45,y=11),label=expression("}"),size=22)+
  geom_text(aes(x=53,y=11),label=expression("16 h"),size=7)+
  
  scale_shape_manual(name="Data from:",values = c("Model" = 18, "Culture" = 1))+
  scale_color_viridis_d(name=expression("PFD"~(μE~"m"^{-2}~"s"^{-1})))+
  theme(
    panel.grid = element_blank(),
    text = element_text(size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.position = c(0.85, 0.45),  # inside plot (x,y)
    legend.background = element_rect(fill = alpha("white", 0.7), colour = NA)
  )
}else{
ggplot()+
  theme_bw()+
  geom_point(data=df_model[df_model$daylight==16&df_model$species=="Alexandrium",],aes(x=Ce*1e9,y=eps,color=(species),shape='Model'),size=4)+
  geom_point(data=df_model[df_model$daylight==24&df_model$species=="Alexandrium",],aes(x=Ce*1e9,y=eps,color=(species),shape='Model'),size=4)+
  geom_smooth(data=df_model[df_model$daylight==16&df_model$species=="Alexandrium",],aes(x=Ce*1e9,y=eps,color=(species)),method='lm',formula=y~poly(x,2),se=F,show.legend = F)+
  geom_smooth(data=df_model[df_model$daylight==24&df_model$species=="Alexandrium",],aes(x=Ce*1e9,y=eps,color=(species)),method='lm',formula=y~poly(x,2),se=F,show.legend = F)+
  geom_point(data=df_model[df_model$PFD==250,],aes(x=Ce*1e9,y=eps,color=(species),shape='Model'),size=4)+
  geom_point(data=df_model[df_model$PFD==55,],aes(x=Ce*1e9,y=eps,color=(species),shape='Model'),size=4)+
  geom_smooth(data=df_model[df_model$PFD==250,],aes(x=Ce*1e9,y=eps,color=(species)),method='lm',formula=y~poly(x,2),se=F,show.legend = F)+
  geom_smooth(data=df_model[df_model$PFD==55,],aes(x=Ce*1e9,y=eps,color=(species)),method='lm',formula=y~poly(x,2),se=F,show.legend = F)+
  xlab(expression(CO[2~"("*aq*")"]*" (μmol/L)"))+
  ylab(expression(epsilon[p]~"(‰)"))+
  
  geom_text(aes(x=70,y=23),label=expression("}"),size=18)+
  geom_text(aes(x=78,y=23),label=expression("24 h"),size=7)+
  geom_text(aes(x=58,y=14),label=expression("}"),size=22)+
  geom_text(aes(x=66,y=14),label=expression("16 h"),size=7)+
  
  
  geom_point(data=culture.data,aes(x=(co2*1.025),y=eps,color=(species),shape='Culture'),size=4)+
  scale_shape_manual(name="Data from:",values = c("Model" = 18, "Culture" = 1))+
  scale_color_manual(name="Species",values = c("red","blue","purple","orange"),labels=c(expression(italic("A. tamarense")),expression(italic("G. spinifera")),expression(italic("P. reticulatum")),expression(italic("S. trochoidea"))))+
  ylim(0,26.7)+
  xlim(0,100)+
  theme(
    panel.grid = element_blank(),
    text = element_text(size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.position = c(0.89, 0.37),  # inside plot (x,y)
    legend.background = element_rect(fill = alpha("white", 0.7), colour = NA)
  )
}

######################################
##Plot relative HCO3- uptake vs pCO2
######################################
if(unique(culture.data$species=="E. huxleyi")){
  
ggplot()+
  geom_point(data=df_model,aes(x=pCO2,y=Uptake,color="Model"),size=5,shape=15,show.legend = F)+
  theme_bw()+
  geom_point(data=fluxdata,aes(x=co2,y=Uptake,color="Culture"),size=5,shape=15,show.legend = F)+
  geom_errorbar(data=fluxdata,aes(x=co2,y=Uptake,ymin=Uptake-Uptakesd,ymax=Uptake+Uptakesd),show.legend = F)+
  xlab(expression(italic(p)*CO[2]~"(µatm)"))+
  ylab(expression("HCO"[3]^"-"~"/total DIC uptake"))+
  scale_shape_discrete(name="Daylength (h)")+
  scale_color_viridis_d(name="Data from:")+
  ylim(c(0,1))+
  theme(text=element_text(size=20),
        panel.grid = element_blank(),
        legend.text = element_text(size=15),legend.title = element_text(size=15))

}else{
  ggplot()+
    geom_point(data=df_model[df_model$PFD>55&df_model$daylight==16,],aes(x=Ce*1e9,y=Uptake,shape=species,color="Model"),size=4,show.legend = F)+
    geom_point(data=fluxdata,aes(x=co2,y=Uptake,color="Culture",shape=Species),size=4,show.legend = F)+
    geom_errorbar(data=fluxdata,aes(x=co2,ymin=Uptake-4.3*Uptakesd/1.73,ymax=Uptake+4.3*Uptakesd/1.73,color="Culture"),width=0,show.legend = F)+
    theme_bw()+
    scale_color_viridis_d(name="")+
    ylim(c(0,1))+
    xlab(expression(CO[2~"("*aq*")"]*" (μmol/L)"))+
    ylab(expression("HCO"[3]^"-"~"/total DIC uptake"))+
    scale_shape_discrete(name="Species",labels=c(expression(italic("A. tamarense")),expression(italic("G. spinifera")),expression(italic("P. reticulatum")),expression(italic("S. trochoidea"))))+
    theme(text=element_text(size=20),
          panel.grid = element_blank(),
          legend.text = element_text(size=15),legend.title = element_text(size=15))
  
}
######################################
##Plot CO2 leakage vs pCO2
######################################
if(unique(culture.data$species=="E. huxleyi")){
  
ggplot()+
  geom_point(data=df_model,aes(x=pCO2,y=Leakage,color="Model"),size=4,shape=15)+
  geom_point(data=fluxdata,aes(x=co2,y=leakage,color="Culture"),size=4,shape=15)+
  theme_bw()+
  xlab(expression(italic(p)*CO[2]~(mu*atm)))+
  ylab(expression("Leakage"))+
  scale_color_viridis_d(name="")+
  ylim(c(0,1))+
  scale_shape_discrete(name="Daylength (h)")+
  theme(text=element_text(size=20),
        panel.grid = element_blank(),
        legend.position = c(0.8,0.8),
        legend.text = element_text(size=15),legend.title = element_text(size=15))
}else{
  ggplot()+
    geom_point(data=df_model[df_model$PFD>55&df_model$daylight==16,],aes(x=Ce*1e9,y=Leakage,shape=species,color="Model"),size=5)+
    geom_point(data=fluxdata,aes(x=co2,y=leakage,color="Culture",shape=Species),size=4)+
    geom_errorbar(data=fluxdata,aes(x=co2,ymin=leakage-4.3*leakagesd/1.73,ymax=leakage+4.3*leakagesd/1.73,color="Culture"),width=0,show.legend = F)+
    theme_bw()+
    xlab(expression(CO[2~"("*aq*")"]*" (μmol/L)"))+
    ylab("Leakage")+
    ylim(c(0,1))+
    scale_color_viridis_d(name="")+
    scale_shape_discrete(name="Species",labels=c(expression(italic("A. tamarense")),expression(italic("G. spinifera")),expression(italic("P. reticulatum")),expression(italic("S. trochoidea"))))+
    theme(text=element_text(size=20),
          panel.grid = element_blank(),
          legend.position = c(0.8,0.25),
          legend.text = element_text(size=15),legend.title = element_text(size=15))
  
}

######################################
##Plot Epsilon_PIC
######################################
if(unique(culture.data$species=="E. huxleyi")){
  
ggplot(data=df_model[df_model$ref=="Rost et al. 2002",],aes(x=PFD,y=EpcalcCO2,color=(Ce*1e9),shape="Model"))+
  geom_point(data=culture.data[culture.data$ref=="Rost et al. 2002",],aes(x=PFD,y=EpcalcCO2,shape="Culture",color=co2*1.025),size=6)+
  geom_point(size=4)+
  theme_bw()+
  scale_shape_manual(name="Data from:",values = c("Model" = 18, "Culture" = 1))+
  xlab(expression(PFD~(μE~m^{-2}~s^{-1})))+
  ylab(expression(epsilon[PIC]~"(‰)"))+
  scale_color_viridis_c(name=expression(CO[2]~(mu*mol~L^{-1})))+
  theme(text=element_text(size=15),panel.grid = element_blank(),
        legend.position = c(0.7,0.6),legend.background = element_blank())
}
