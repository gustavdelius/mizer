# Simulations. Exploring the community effects of anchovy's fishing mortality
# ----------------------------------------------
# 3) Running the simulations  --- project()
# ----------------------------------------------
# Scoping the background information for an ecosystem approach to fisheries in Scottish waters:
# Review of predator-prey interactions with fisheries, and balanced harvesting


rm(list = ls())
library(mizer)
setwd("~/WORK_MLLA/2018_POSTDOC_PUC/Objetivo3_MIZER_NCME/01_NCME_MIZER/NCME_MIZER_V2/mizer/inst/humboldt/")

# Second set of parameters for an steady state fish community 
# -----------------------------------------------------------
# habia un *13.rds file
humboldt_sscomm <- readRDS("~/WORK_MLLA/2018_POSTDOC_PUC/Objetivo3_MIZER_NCME/01_NCME_MIZER/NCME_MIZER_V2/mizer/inst/humboldt/params (7)_v12.rds")
sscom_inicon_sim=humboldt_sscomm

# summary(sscom_inicon_sim)
# An object of class "MizerParams" 
# Community size spectrum:
#   minimum size:	1e-04
# maximum size:	573889
# no. size bins:	400
# Background size spectrum:
#   minimum size:	1.02453e-09
# maximum size:	573889
# no. size bins:	60
# Species details:
#   species    w_inf     w_mat   beta sigma
# Mesopelagic Mesopelagic      2.2      0.07 105.00  2.05
# Anchovy         Anchovy     66.5     14.41   0.33  5.36
# Sardine         Sardine    625.1    181.51 205.00  5.08
# Mackerel       Mackerel   2008.3    353.54 163.00  3.07
# JMackerel     JMackerel   4554.0    200.56 144.00  2.39
# Palm_Ruff     Palm_Ruff  12178.0   1494.43  99.00  1.14
# EP_Bonito     EP_Bonito  13333.1   1693.34  33.00  1.27
# Swordfish     Swordfish 573888.7 132683.60 652.00  1.65
# Fishing gear details:
#   Gear			Target species
# sigmoid_gear 		 Sardine Mackerel JMackerel Palm_Ruff EP_Bonito Swordfish 
# sigmoid_gear_Anchovy 		 Anchovy 

# ---------------------------------------------
# Initial Community at Steady states
# ---------------------------------------------
sim0<- project(sscom_inicon_sim, effort = 1, t_max =50,  t_save = 1)

sim0@params@linecolour["Plankton"]    <- "darkgreen"
sim0@params@linecolour["Mesopelagic"] <- "orange"
sim0@params@linecolour["Anchovy"]     <- 'blue'
sim0@params@linecolour["Sardine"]     <- 'red'
sim0@params@linecolour["Mackerel"]    <- 'limegreen'
sim0@params@linecolour["JMackerel"]   <- 'orange'
sim0@params@linecolour["Palm_Ruff"]   <- 'blue'
sim0@params@linecolour["EP_Bonito"]   <- 'red'
sim0@params@linecolour["Swordfish"]   <- 'limegreen'

sim0@params@linetype["Plankton"]    <- "solid"
sim0@params@linetype["Mesopelagic"] <- "solid"
sim0@params@linetype["Anchovy"] <- "solid"
sim0@params@linetype["Sardine"] <- "solid"
sim0@params@linetype["Mackerel"] <- 'solid'
sim0@params@linetype["JMackerel"] <- 'dashed'
sim0@params@linetype["Palm_Ruff"] <- 'dashed'
sim0@params@linetype["EP_Bonito"] <- 'dashed'
sim0@params@linetype["Swordfish"] <- 'dashed'

plot(sim0)

# ---- Species information ----
Nf_ss0=sim0@n
Bf_ss0=getBiomass(sim0)
Np_ss0=sim0@n_pp
M_ss0=getM2(sim0)  
N_sp0=getN(sim0)
BS_sp0=getSSB(sim0)
Y_sp0=getYield(sim0)

# ---- Community information ----
CommSp0=getCommunitySlope(sim0)
CommMMW0=getMeanMaxWeight(sim0)
CommMeanW0=getMeanWeight(sim0)

output_S0<-list(Nf_ss0,Bf_ss0,Np_ss0,M_ss0,N_sp0,BS_sp0,Y_sp0,CommSp0,CommMMW0,CommMeanW0)

save(output_S0,file='SIM_report_S0.Rdata')
load("SIM_report_S0.Rdata")      

# ---------------------------------------------------------------------------
# 2. Setting constant effort for different gears (E rescaled the selectivity)
# ---------------------------------------------------------------------------
# Scenario Fanc goes up

Effort <- c(sigmoid_gear = 1, sigmoid_gear_Anchovy = 0.5)
sim <- project(sscom_inicon_sim, effort = Effort, t_max = 50, dt = 0.1, t_save = 1)
plot(sim)

# ---- Species information ----
Nf_ss=sim@n
Bf_ss=getBiomass(sim)
Np_ss=sim@n_pp
M_ss=getM2(sim)  
N_sp=getN(sim)
BS_sp=getSSB(sim)
Y_sp=getYield(sim)

# ---- Community information ----
CommSp=getCommunitySlope(sim)
CommMMW=getMeanMaxWeight(sim)
CommMeanW=getMeanWeight(sim)

output_S2<-list(Nf_ss,Bf_ss,Np_ss,M_ss,N_sp,BS_sp,Y_sp,CommSp,CommMMW,CommMeanW)

save(output_S2,file='SIM_report_S2.Rdata')
load("SIM_report_S2.Rdata")      

# Results plots - M Canales Report Fondecyt
# ----------------------------------------------

# Validation model. Two graphs Growth curve and catch size distribution

# Plots of the Growth Curve
# --------------------------
# Read VB parameters first

rm(list = c())
library(mizer)
library(plotly)
packageVersion('plotly')
library(ggplot2)

# Second set of parameters for an steady state fish community 
# -----------------------------------------------------------
humboldt_sscomm <- readRDS("~/WORK_MLLA/2018_POSTDOC_PUC/Objetivo3_MIZER_NCME/01_NCME_MIZER/NCME_MIZER_V2/mizer/inst/humboldt/params (7)_v12.rds")

m_growthC <- getGrowthCurves(humboldt_sscomm)
parmod_Mlla <- read.csv("/Users/Mlla/WORK_MLLA/2018_POSTDOC_PUC/Objetivo3_MIZER_NCME/01_NCME_MIZER/NCME_MIZER_V2/mizer/inst/humboldt/speciesNCME_Mariella.csv")
Winf_vb<-parmod_Mlla$w_inf[c(2,4,6,8,10,12,14,16)]
Linf_vb<-parmod_Mlla$Linf[c(2,4,6,8,10,12,14,16)]
k_vb<-parmod_Mlla$k_vb[c(2,4,6,8,10,12,14,16)]
t_vb<-parmod_Mlla$to[c(2,4,6,8,10,12,14,16)]
b2<-parmod_Mlla$b2[c(2,4,6,8,10,12,14,16)]

com.age=as.numeric(colnames(m_growthC))
vb_gc=array(0,dim=c(8,length(com.age)))

for (i in 1:8){
  for (j in 1:length(com.age))
  {
    vb_gc[i,j]=Winf_vb[i]*(1-exp(-k_vb[i]*(com.age[j]-t_vb[i])))^2
  }
}

m_gC<-as.data.frame(m_growthC)
vb_gC<-as.data.frame(vb_gc)

rownames(vb_gC)<-c("Mesopelagic","Anchovy","Sardine","Mackerel",'JMackerel','Palm_Ruff', 'EP_Bonito', "Swordfish")
colnames(vb_gC)<- as.character(com.age)
g1=melt(t(m_gC)) 
g1$gr_type<-'Model'
g2=melt(t(vb_gC))
g2$gr_type<-'Observed'

dat=rbind(g1,g2)
colnames(dat)=c('age','sps','Growth','Legend')

p <- ggplot(dat, aes(x=age, y = Growth),size=1) +
  geom_line(aes(linetype=Legend)) +
  scale_linetype_manual(values=c("solid",'dashed')) +
  facet_wrap(~sps, scales = "free", nrow = 2) + theme_bw() +
  labs(x = "Age (years)")+
  labs(y = "Somatic Growth (g)") +
  theme(axis.title = element_text( size = 14, face = "bold" ),
        legend.position="top",
        # The new stuff
        strip.text = element_text(size = 15),
        axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        axis.title.x = element_text(size=18,face='bold'),
        axis.title.y = element_text(size=18,face='bold')) 
p

# Plots Catch size distributions
# ------------------------------

# AQUI PLOTEAR ESTA VALADACION








# Individual species spawning biomass plots by scenarios
# ------------------------------------------------------
rm(list=ls())

load("SIM_report_S0.Rdata")
load("SIM_report_S1.Rdata")
load("SIM_report_S2.Rdata")
load("SIM_report_S3.Rdata")
load("SIM_report_S4.Rdata")

# --- 
sb_S0=as.data.frame(output_S0[[6]])
sb_S1=as.data.frame(output_S1[[6]])
sb_S2=as.data.frame(output_S2[[6]])
sb_S3=as.data.frame(output_S3[[6]])
sb_S4=as.data.frame(output_S4[[6]])

# --- 
years <- 0:50
I1_S0 <-cbind(years,(sb_S0/sb_S0))
I1_S1 <-cbind(years,(sb_S1/sb_S0))
I1_S2 <-cbind(years,(sb_S2/sb_S0))
I1_S3 <-cbind(years,(sb_S3/sb_S0))
I1_S4 <-cbind(years,(sb_S4/sb_S0))

library(reshape2)
sb0_df <- melt(I1_S0,id.vars=c("years"))
sb1_df <- melt(I1_S1,id.vars=c("years"))
sb2_df <- melt(I1_S2,id.vars=c("years"))
sb3_df <- melt(I1_S3,id.vars=c("years"))
sb4_df <- melt(I1_S4,id.vars=c("years"))

sb_df <- rbind(
  cbind(sb0_df, scenario = "S0"),
  cbind(sb1_df, scenario = "S1"),
  cbind(sb2_df, scenario = "S2"),
  cbind(sb3_df, scenario = "S3"),
  cbind(sb4_df, scenario = "S4"))

colnames(sb_df)<-c("time","sp",'Spawning_Biomass','Scenarios')

library(ggplot2)
p <- ggplot(sb_df, aes(x=time, y =Spawning_Biomass)) + 
  geom_line(aes(color=Scenarios,linetype=Scenarios),size=1) +
  scale_linetype_manual(values=c("solid",'solid','dotted','solid','dotted')) +
  scale_color_manual(values=c("red",'black', 'black','green','green')) +
  facet_wrap(~sp, scales = "free", nrow = 2) + theme_bw() +
  theme(axis.title = element_text( size = 14, face = "bold"),
        legend.position="top",
        strip.text = element_text(size = 15),
        axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        axis.title.x = element_text(size=14,face='bold'),
        axis.title.y = element_text(size=14,face='bold')) +
  labs(x="Years", y=expression("Relative spawning biomass [gm"^-3~']')) 
p


# Individual species Yield plots and scenarios
# ----------------------------------------------
rm(list=ls())

load("SIM_report_S0.Rdata")
load("SIM_report_S1.Rdata")
load("SIM_report_S2.Rdata")
load("SIM_report_S3.Rdata")
load("SIM_report_S4.Rdata")

# --- 
Y_S0=as.data.frame(output_S0[[7]])
Y_S1=as.data.frame(output_S1[[7]])
Y_S2=as.data.frame(output_S2[[7]])
Y_S3=as.data.frame(output_S3[[7]])
Y_S4=as.data.frame(output_S4[[7]])

# --- 
years <- 0:50
I1_S0 <-cbind(years,(Y_S0/Y_S0))
I1_S0$Mesopelagic<-0 
I1_S1 <-cbind(years,(Y_S1/Y_S0))
I1_S1$Mesopelagic<-0 
I1_S2 <-cbind(years,(Y_S2/Y_S0))
I1_S2$Mesopelagic<-0 
I1_S3 <-cbind(years,(Y_S3/Y_S0))
I1_S3$Mesopelagic<-0 
I1_S4 <-cbind(years,(Y_S4/Y_S0))
I1_S4$Mesopelagic<-0 

library(reshape2)
yb0_df <- melt(I1_S0,id.vars=c("years"))
yb1_df <- melt(I1_S1,id.vars=c("years"))
yb2_df <- melt(I1_S2,id.vars=c("years"))
yb3_df <- melt(I1_S3,id.vars=c("years"))
yb4_df <- melt(I1_S4,id.vars=c("years"))

yb_df <- rbind(
  cbind(yb0_df, scenario = "S0"),
  cbind(yb1_df, scenario = "S1"),
  cbind(yb2_df, scenario = "S2"),
  cbind(yb3_df, scenario = "S3"),
  cbind(yb4_df, scenario = "S4"))

colnames(yb_df)<-c("time","sp",'Yield','Scenarios')

library(ggplot2)
p <- ggplot(yb_df, aes(x=time, y = Yield)) + 
  geom_line(aes(color=Scenarios,linetype=Scenarios),size=1) +
  scale_linetype_manual(values=c("solid",'solid','dotted','solid','dotted')) +
  scale_color_manual(values=c("red",'black', 'black','green','green')) +
  facet_wrap(~sp, scales = "free", nrow = 2) + theme_bw() +
  theme(axis.title = element_text( size = 14, face = "bold"),
        legend.position="top",
        strip.text = element_text(size = 15),
        axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        axis.title.x = element_text(size=16,face='bold'),
        axis.title.y = element_text(size=16,face='bold')) +
  labs(x="Years", y=expression("Relative yield [gm"^-3~']')) 
p

# Community plots by scenarios
# ----------------------------------------------
rm(list = ls())
library(mizer)
setwd("~/WORK_MLLA/2018_POSTDOC_PUC/Objetivo3_MIZER_NCME/01_NCME_MIZER/NCME_MIZER_V2/mizer/inst/humboldt/")

load("SIM_report_S0.Rdata")
load("SIM_report_S1.Rdata")
load("SIM_report_S2.Rdata")
load("SIM_report_S3.Rdata")
load("SIM_report_S4.Rdata")

# --- 
I_S0=cbind(as.data.frame(output_S0[[9]])[2],as.data.frame(output_S0[[10]]),as.data.frame(output_S0[[8]][,1]))
colnames(I_S0)<-c('wmm',"wm","slope")
I_S1=cbind(as.data.frame(output_S1[[9]])[2],as.data.frame(output_S1[[10]]),as.data.frame(output_S1[[8]][,1]))
colnames(I_S1)<-c('wmm',"wm","slope")
I_S2=cbind(as.data.frame(output_S2[[9]])[2],as.data.frame(output_S2[[10]]),as.data.frame(output_S2[[8]][,1]))
colnames(I_S2)<-c('wmm',"wm","slope")
I_S3=cbind(as.data.frame(output_S3[[9]])[2],as.data.frame(output_S3[[10]]),as.data.frame(output_S3[[8]][,1]))
colnames(I_S3)<-c('wmm',"wm","slope")
I_S4=cbind(as.data.frame(output_S4[[9]])[2],as.data.frame(output_S4[[10]]),as.data.frame(output_S4[[8]][,1]))
colnames(I_S4)<-c('wmm',"wm","slope")

# --- 
years <- 0:50
wmm0_S0 <-as.data.frame(cbind(years,(I_S0[,1]/I_S0[,1]),(I_S0[,2]/I_S0[,2]),(I_S0[,3]/I_S0[,3])))
colnames(wmm0_S0)<-c('years','Mean maximun weight','Mean weight','Slope')

wmm1_S1 <-as.data.frame(cbind(years,(I_S1[,1]/I_S0[,1]),(I_S1[,2]/I_S0[,2]),(I_S1[,3]/I_S0[,3])))
colnames(wmm1_S1)<-c('years','Mean maximun weight','Mean weight','Slope')

wmm2_S2 <-as.data.frame(cbind(years,(I_S2[,1]/I_S0[,1]),(I_S2[,2]/I_S0[,2]),(I_S2[,3]/I_S0[,3])))
colnames(wmm2_S2)<-c('years','Mean maximun weight','Mean weight','Slope')

wmm3_S3 <-as.data.frame(cbind(years,(I_S3[,1]/I_S0[,1]),(I_S3[,2]/I_S0[,2]),(I_S3[,3]/I_S0[,3])))
colnames(wmm3_S3)<-c('years','Mean maximun weight','Mean weight','Slope')

wmm4_S4 <-as.data.frame(cbind(years,(I_S4[,1]/I_S0[,1]),(I_S4[,2]/I_S0[,2]),(I_S4[,3]/I_S0[,3])))
colnames(wmm4_S4)<-c('years','Mean maximun weight','Mean weight','Slope')

library(reshape2)
wmm0_S0 <- melt(wmm0_S0,id.vars=c("years"))
wmm1_S1 <- melt(wmm1_S1,id.vars=c("years"))
wmm2_S2 <- melt(wmm2_S2,id.vars=c("years"))
wmm3_S3 <- melt(wmm3_S3,id.vars=c("years"))
wmm4_S4 <- melt(wmm4_S4,id.vars=c("years"))

wmm_df <- rbind(
  cbind(wmm0_S0, scenario = "S0"),
  cbind(wmm1_S1, scenario = "S1"),
  cbind(wmm2_S2, scenario = "S2"),
  cbind(wmm3_S3, scenario = "S3"),
  cbind(wmm4_S4, scenario = "S4"))

colnames(wmm_df)<-c("Year","type_w",'Values','Scenarios')
#ubi<-which(wmm_df$type_w!='Slope')
#wmm_df<-wmm_df[ubi,]

library(ggplot2)
p <- ggplot(wmm_df, aes(x=Year, y =Values)) + 
  geom_line(aes(color=Scenarios,linetype=Scenarios),size=1) +
  scale_linetype_manual(values=c("solid",'solid','dotted','solid','dotted')) +
  scale_color_manual(values=c("red",'black', 'black','green','green')) +
  facet_wrap(~type_w, scales = "free", nrow = 1) + theme_bw() +
  theme(axis.title = element_text( size = 14, face = "bold"),
        legend.position="top",
        strip.text = element_text(size = 15),
        axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        axis.title.x = element_text(size=16,face='bold'),
        axis.title.y = element_text(size=16,face='bold')) +
  labs(x="Years", y=expression("Relative weight [g]")) 
p


# Individual species M2 plots
# ---------------------------
# This is telling me how much of each predator is consuming at a particular body size
# predRate_sps<-(getPredRate(humboldt_sscomm,humboldt_sscomm@initial_n, humboldt_sscomm@initial_n_pp))
# bodysize_com<-humboldt_sscomm@dw_full
# bodysize_fish<-humboldt_sscomm@dw


