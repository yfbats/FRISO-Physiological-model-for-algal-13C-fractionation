#This script can be used to investigate generic model behavior based on generic artificial data
#stored to 'yourdata'

###################################
##Set correct working directory
###################################
#setwd("")

######################################
##Required packages
######################################
packages <- c("seacarb","readxl","ggplot2","dplyr")
lapply(packages, require, character.only = TRUE)

################################################################################
##Create artificial data 
##In this example CO2 and PFD are varying, and the other parameters are constant
################################################################################
co2_vec=rep(seq(1,40,by=1),3)
PFD_vec=rep(c(50,100,150),each=40)

yourdata = data.frame(co2_vec,                                        #μM
                    PFD_vec,                                          #μE/m2/s
                    
                    pH=rep(8,length(co2_vec*PFD_vec)),
                    Vcell=rep(1e-10,length(co2_vec*PFD_vec)),         #cm3
                    Temp=rep(18,length(co2_vec*PFD_vec)),             #Celsius
                    sal=rep(35,length(co2_vec*PFD_vec)),              #PSU
                    r1=rep(1e-8,length(co2_vec*PFD_vec)),             #mol/cm3/s
                    daylight=rep(16,length(co2_vec*PFD_vec)),
                    species=rep("E. huxleyi",length(co2_vec*PFD_vec)) #Species settings
                  )

########################################
##run model and save results in df_model
########################################
source("Carbon flux and isotope model.R")
df_model=solver_model(yourdata)


###################################
##Plot model output
###################################

ggplot()+
  theme_bw()+
  ylim(0,26.7)+
  xlim(0,90)+
  geom_point(data=df_model,aes(x=Ce*1e9,y=eps,color=as.factor(PFD),shape='Model'),size=4)+
  xlab(expression(CO[2~"("*aq*")"]*" (μmol/L)"))+
  ylab(expression(epsilon[p]~"(‰)"))+
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
