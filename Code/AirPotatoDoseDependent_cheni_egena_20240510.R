library(tidyverse)
library(car)
library(emmeans)
library(glmmTMB)
library(patchwork)
library(performance)
library(viridis)
library(DHARMa)
library(GGally)
library(skimr)
# # Load cheni dose dependent data file -----------------------------------

ddc<-read.csv("Data/Cheni_dose_dependent_assay.csv")


# Explore and process the file --------------------------------------------

ddc1 <- ddc %>% mutate(Beetle.Wgt=Beetle.FW-Beetle.IW,
                    Frass.Wgt=Frass.tube.FW-Frass.tube.IW,
                    Treatment=fct_relevel(Treatment, c("Control","2.5mg","5mg","10mg","15mg","20mg")),
                    Trt=case_when(Treatment == 'Control' ~ 0,
                                  Treatment =='2.5mg' ~ 2.5,
                                  Treatment == '5mg' ~ 5.0,
                                  Treatment == '10mg' ~ 10,
                                  Treatment == '15mg' ~ 15,
                                  Treatment == '20mg' ~ 20),
                    Trt=as.numeric(Trt), Surv=(1-Dead)) %>% 
  filter(Frass.Wgt<.1) %>% #remove 1 very high frass weight and it also died early and had no final weight of beetle
  drop_na(Beetle.Wgt)  #all the ones that were dead and few missing

head(ddc1) 



ddc2 <- ddc %>% mutate(Beetle.Wgt=Beetle.FW-Beetle.IW,
                          Frass.Wgt=Frass.tube.FW-Frass.tube.IW,
                          Treatment=fct_relevel(Treatment, c("Control","2.5mg","5mg","10mg","15mg","20mg")),
                          Trt=case_when(Treatment == 'Control' ~ 0,
                                        Treatment =='2.5mg' ~ 2.5,
                                        Treatment == '5mg' ~ 5.0,
                                        Treatment == '10mg' ~ 10,
                                        Treatment == '15mg' ~ 15,
                                        Treatment == '20mg' ~ 20),
                          Trt=as.numeric(Trt), Surv=(1-Dead)) %>% 
                          filter(Frass.Wgt<.1)

survival_summary_cheni <- ddc2 %>%
  group_by(Treatment) %>%
  summarise(
    Total = n(),
    Survived = sum(Surv),
    Died = Total - Survived)


ddc_resp <- ddc1 %>% select(Trt, Beetle.Wgt, Frass.Wgt, Surv )


summary(ddc_resp)
skim(ddc_resp)
ggpairs(ddc_resp) #frass and beetle weight are moderately positively correlated, frass and survival are strong positive correlated


ch1 <- lm(cbind(Beetle.Wgt, Frass.Wgt, Surv)~ Trt, data= ddc1)
Manova(ch1, test.statistic= "Wilks") #P=0.1319 not significant

# Testing different models for each individual variable!!! ####
# Beetle weight -----------------------------------------------------------
lm1ch1 <- glmmTMB(Beetle.Wgt~1, data=ddc1) # null model
lm1ch2 <- glmmTMB(Beetle.Wgt~Trt, data=ddc1) # linear model
lm1ch3 <- glmmTMB(Beetle.Wgt~Trt+I(Trt^2), data=ddc1) # quadratic model 
anova(lm1ch1,lm1ch2,lm1ch3) # quadratic model anova p =  0.0901 . lower aic than linear model 694.28

Anova(lm1ch3) # trt p = 0.07530 . ;  trt^2 p = 0.08804 .
summary(lm1ch3)
simulateResiduals(lm1ch3, plot=T)
MuMIn::r.squaredGLMM(lm1ch3) # R2m= 0.02772571 ; R2c= 0.02772571

ddc1ae <- ggplot(ddc1%>% filter(Beetle.Wgt<20), aes(x=Trt, y=Beetle.Wgt, fill=Trt))+
  geom_smooth(method="lm", formula = y~x+I(x^2), color="#395D9C", fill="#395D9C", lwd=2)+
  geom_jitter(height=0, width=.1, pch=21, size=3) +
  scale_fill_viridis(direction=-1, option="G")+
  labs(y="Change in beetles mass (mg)", x="Diosgenin in diet (mg/g)")+
  theme_bw(base_size = 24)+
  theme(legend.position = "none"); ddc1ae

# Frass weight -----------------------------------------------------------
lm2ch1 <- glmmTMB(Frass.Wgt*1000~1, data=ddc1)
lm2ch2 <- glmmTMB(Frass.Wgt*1000~Trt, data=ddc1)
lm2ch3 <- glmmTMB(Frass.Wgt*1000~Trt+I(Trt^2), data=ddc1)
anova(lm2ch1,lm2ch2,lm2ch3) 
# with NAs quadratic model anova p = 0.001148 ** lower aic than linear model -869.48
# without NAs quadratic model anova p = 0.05977 ** lower aic than linear model -719.03

Anova(lm2ch3) # trt p=0.0004066 ***;  trt^2 p= 0.0009242 *** || without NAs trt^2 p=0.0578;  trt p= 0.0174 
summary(lm2ch3)
simulateResiduals(lm2ch3, plot=T)
MuMIn::r.squaredGLMM(lm2ch3) 
# WIth NAs R2m, R2c= 0.08174868 
#Without NAs R2m, R2c= 0.06693664 


ddc2be <- ggplot(ddc1 , aes(x=Trt, y=Frass.Wgt*1000, fill=Trt))+
  geom_smooth(method="lm", formula = y~x+I(x^2), color="#395D9C", fill="#395D9C", lwd=2)+
  geom_jitter(height=0, width=.1, pch=21, size=3) +
  scale_fill_viridis(direction=-1, option="G")+
  labs(y="Frass mass (mg)", x="Diosgenin in diet (mg/g)")+
  theme_bw(base_size = 24)+
  theme(legend.position = "none"); ddc2be


# Survival -----------------------------------------------------------
lm3ch1 <- glmmTMB(Surv~1, data=ddc2, family=binomial )
lm3ch2 <- glmmTMB(Surv~Trt, data=ddc2, family=binomial )
lm3ch3 <- glmmTMB(Surv~Trt+I(Trt^2), data=ddc2, family=binomial )
anova(lm3ch1,lm3ch2,lm3ch3) # quadratic model anova p = 0.003377 ** lower aic than linear model 141.48

Anova(lm3ch3)  # trt p=0.005010 ;  trt^2 p= 0.005081 **
summary(lm3ch3)
simulateResiduals(lm3ch3, plot=T)
MuMIn::r.squaredGLMM(lm3ch3) #theoretical R2m, R2c= 0.10891973 ; delta R2m, R2c= 0.06104416

ddc2c <- ggplot(ddc2, aes(x=Trt, y=Surv, fill=Trt))+
  geom_smooth(method="glm", formula = y~x+I(x^2), method.args=list(family=binomial),
              color="#395D9C", fill="#395D9C", lwd=2)+
  geom_jitter(height=0.05, width=.25, pch=21, size=2) +
  scale_fill_viridis(direction=-1, option="G")+
  scale_y_continuous(labels = scales::percent)+
  labs(y="Survival", x="Diosgenin in diet (mg/g)")+
  theme_bw(base_size = 24)+
  theme(legend.position = "none"); ddc2c


beetle_plot_cheni_dose <- (ddc2be+ddc2c) +
  plot_annotation(tag_levels = 'a') + 
  plot_layout(guides='collect') & theme(legend.position='none');beetle_plot_cheni_dose

ggsave("beetle_plot_cheni_dose2.tiff", beetle_plot_cheni_dose, width=12, height=6, units="in", dpi=600, compression = "lzw", path="Outputs")



# # Load egena dose dependent data file -----------------------------------
dd<-read.csv("Data/Egena_dose_dependent_assay.csv") 
head(dd)  

# Explore and process the file --------------------------------------------

unique(dd$Treatment)

dd <- dd %>% mutate(Beetle.Wgt=Beetle.FW-Beetle.IW,
                      Frass.Wgt=Frass.tube.FW-Frass.tube.IW,
                    Trt = fct_relevel(as.factor(Treatment), c("0", "0.15", "0.25", "0.35", "0.5")),
                      Trt=case_when(Treatment == '0' ~ 0,
                                    Treatment =='0.15' ~ 0.15,
                                    Treatment == '0.25' ~ 0.25,
                                    Treatment == '0.35' ~ 0.35,
                                    Treatment == '0.5' ~ 0.5),
                      Trt=as.numeric(Trt), Surv=(1-Dead))

head(dd) 

survival_summary_egena <- dd %>%
  group_by(Treatment) %>%
  summarise(
    Total = n(),
    Survived = sum(Surv),
    Died = Total - Survived)


dd_resp <- dd %>% select(Trt, Beetle.Wgt, Frass.Wgt, Surv )
summary(dd_resp)
skim(dd_resp)
ggpairs(dd_resp) 
#trt-frass, trt-weight, trt- surv -- all negatively related, 
# beetle and frass- mod positive
# beetle and surv-- weak positive
#frass-surv--mod positive

eg1 <- lm(cbind(Beetle.Wgt, Frass.Wgt, Surv)~ Trt, data= dd)
Manova(eg1, test.statistic= "Wilks") #P=0.04376 *

# Testing different models for each variable individually!!! ####

# Beetle weight -----------------------------------------------------------

lm1a1 <- glmmTMB(Beetle.Wgt~1, data=dd) # null model
lm1a2 <- glmmTMB(Beetle.Wgt~Treatment, data=dd) # linear model
lm1a3 <- glmmTMB(Beetle.Wgt~Treatment+I(Treatment^2), data=dd) # quadratic model 
anova(lm1a1,lm1a2,lm1a3) #no p value significant but linear model has lower aic= 1141.3 and lower P = 0.2712

Anova(lm1a2) #P = 0.2712
summary(lm1a2) #P = 0.2712, estimate -16.22, SE 14.704, z = -1.103; intercept estimate 26.645
simulateResiduals(lm1a2, plot=T)
MuMIn::r.squaredGLMM(lm1a2) #R2m, R2c= 0.01012177

dd1ae <- ggplot(dd%>% filter(Beetle.Wgt<20), aes(x=Treatment*100, y=Beetle.Wgt*100, fill=Treatment))+
  geom_smooth(method="lm", formula = y~x, color="#395D9C", fill="#395D9C", lwd=2)+
  geom_jitter(height=0, width=.1, pch=21, size=3) +
  scale_fill_viridis(direction=-1, option="G")+
  labs(y="Change in beetles mass (mg)", x="Diosgenin in diet (mg/g)")+
  theme_bw(base_size = 24) + 
  theme(legend.position = "none"); dd1ae

# Frass weight -----------------------------------------------------------

lm2a1 <- glmmTMB(Frass.Wgt~1, data=dd)
lm2a2 <- glmmTMB(Frass.Wgt~Treatment, data=dd)
lm2a3 <- glmmTMB(Frass.Wgt~Treatment+I(Treatment^2), data=dd)
anova(lm2a1,lm2a2,lm2a3) # linear model has lower aic= -694.06 and lower P = 0.01733 *
simulateResiduals(lm2a2, plot=T)
Anova(lm2a2) #P= 0.01604 *
summary(lm2a2)
MuMIn::r.squaredGLMM(lm2a2) #R2m, R2c= 0.0464632

dd2be <- ggplot(dd , aes(x=Treatment*100, y=Frass.Wgt*1000, fill=Treatment))+
  geom_smooth(method="lm", formula = y~x, color="#395D9C", fill="#395D9C", lwd=2)+
  geom_jitter(height=0, width=.1, pch=21, size=3) +
  scale_fill_viridis(direction=-1, option="G")+
  labs(y="Frass mass (mg)", x="Diosgenin in diet (mg/g)")+
  theme_bw(base_size = 24)+
  theme(legend.position = "none"); dd2be

# Survival -----------------------------------------------------------
lm3a1 <- glmmTMB(Surv~1, data=dd, family=binomial )
lm3a2 <- glmmTMB(Surv~Trt, data=dd, family=binomial )
lm3a3 <- glmmTMB(Surv~Trt+I(Trt^2), data=dd, family=binomial )
anova(lm3a1,lm3a2,lm3a3) # linear model anova p = 0.01743 lower aic than linear model 76.368
simulateResiduals(lm3a2, plot=T)
Anova(lm3a2)  # trt p= 0.02666
summary(lm3a2)
MuMIn::r.squaredGLMM(lm3a2) # Theoretical R2m, R2c= 0.15622464,  delta R2m, R2c= 0.05197148 

dd2a <- ggplot(dd, aes(x=Trt*100, y=Surv, fill=Trt))+
  geom_smooth(method="glm", formula = y~x, method.args=list(family=binomial),
              color="#395D9C", fill="#395D9C", lwd=2)+
  geom_jitter(height=0.05, width=.25, pch=21, size=2) +
  scale_fill_viridis(direction=-1, option="G")+
  scale_y_continuous(labels = scales::percent)+
  labs(y="Survival", x="Diosgenin in diet (mg/g)")+
  theme_bw(base_size = 24)+
  theme(legend.position = "none"); dd2a

beetle_plot_egena_dose <- (dd2be+dd2a) +
  plot_annotation(tag_levels = 'a') + 
  plot_layout(guides='collect') & theme(legend.position='none');beetle_plot_egena_dose

ggsave("beetle_plot_egena_dose2.tiff", beetle_plot_egena_dose, width=12, height=6, units="in", dpi=600, compression = "lzw", path="Outputs")

