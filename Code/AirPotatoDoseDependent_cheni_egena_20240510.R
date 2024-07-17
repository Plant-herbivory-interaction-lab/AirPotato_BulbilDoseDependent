library(tidyverse)
library(car)
library(emmeans)
library(glmmTMB)
library(patchwork)
library(performance)
library(viridis)

# DATA FOR CHENI ####
ddc<-read.csv("Diosgenin Dose dependent feeding assay.csv") 
head(ddc)  

ddc <- ddc %>% mutate(Beetle.Wgt=Beetle.FW-Beetle.IW,
                    Frass.Wgt=Frass.tube.FW-Frass.tube.IW,
                    Treatment=fct_relevel(Treatment, c("Control","2.5mg","5mg","10mg","15mg","20mg")),
                    Trt=case_when(Treatment == 'Control' ~ 0,
                                  Treatment =='2.5mg' ~ 2.5,
                                  Treatment == '5mg' ~ 5.0,
                                  Treatment == '10mg' ~ 10,
                                  Treatment == '15mg' ~ 15,
                                  Treatment == '20mg' ~ 20),
                    Trt=as.numeric(Trt), Surv=(1-Dead)) %>% 
  filter(Frass.Wgt<.1) #remove 1 very high frass weight

# CONTINUOUS ANALYSIS -- USE THIS !!! ####
lm1ch1 <- glmmTMB(Beetle.Wgt~1, data=ddc) # null model
lm1ch2 <- glmmTMB(Beetle.Wgt~Trt, data=ddc) # linear model
lm1ch3 <- glmmTMB(Beetle.Wgt~Trt+I(Trt^2), data=ddc) # quadratic model 
anova(lm1ch1,lm1ch2,lm1ch3)

Anova(lm1ch3)
summary(lm1ch3)

ddc1ae <- ggplot(ddc%>% filter(Beetle.Wgt<20), aes(x=Trt, y=Beetle.Wgt))+
  geom_smooth(method="lm", formula = y~x+I(x^2))+
  geom_point() +
  labs(y="Change in beetles mass (mg)", x="Diosgenin in diet (mg/g)")+
  theme_bw(base_size = 24)

lm2ch1 <- glmmTMB(Frass.Wgt~1, data=ddc)
lm2ch2 <- glmmTMB(Frass.Wgt~Trt, data=ddc)
lm2ch3 <- glmmTMB(Frass.Wgt~Trt+I(Trt^2), data=ddc)
anova(lm2ch1,lm2ch2,lm2ch3)

Anova(lm2ch3)
MuMIn::r.squaredGLMM(lm2ch3)

ddc2be <- ggplot(ddc , aes(x=Trt, y=Frass.Wgt*1000, fill=Trt))+
  geom_smooth(method="lm", formula = y~x+I(x^2), color="#395D9C", fill="#395D9C", lwd=2)+
  geom_jitter(height=0, width=.1, pch=21, size=3) +
  scale_fill_viridis(direction=-1, option="G")+
  labs(y="Frass mass (mg)", x="Diosgenin in diet (mg/g)")+
  theme_bw(base_size = 24)+
  theme(legend.position = "none")

lm3ch1 <- glmmTMB(Surv~1, data=ddc, family=binomial )
lm3ch2 <- glmmTMB(Surv~Trt, data=ddc, family=binomial )
lm3ch3 <- glmmTMB(Surv~Trt+I(Trt^2), data=ddc, family=binomial )
anova(lm3ch1,lm3ch2,lm3ch3)

Anova(lm3ch3)
MuMIn::r.squaredGLMM(lm3ch3)

ddc2c <- ggplot(ddc, aes(x=Trt, y=Surv, fill=Trt))+
  geom_smooth(method="glm", formula = y~x+I(x^2), method.args=list(family=binomial),
              color="#395D9C", fill="#395D9C", lwd=2)+
  geom_jitter(height=0.05, width=.25, pch=21, size=2) +
  scale_fill_viridis(direction=-1, option="G")+
  scale_y_continuous(labels = scales::percent)+
  labs(y="Survival", x="Diosgenin in diet (mg/g)")+
  theme_bw(base_size = 24)+
  theme(legend.position = "none")

ddcplot2 <- ddc2be+ddc2c
ggsave("ddcPlot2cheni.tiff", ddcplot2, width=6, height=4, units="in", dpi=600, compression = "lzw")


# DATA FOR EGENA ####
dd<-read.csv("Egena_dose_dependent_assay_copy.csv") 
head(dd)  

dd <- dd %>% mutate(Beetle.Wgt=Beetle.FW-Beetle.IW,
                    Frass.Wgt=Frass.tube.FW-Frass.tube.IW)


# CONTINUOUS ANALYSIS -- USE THIS !!! ####
lm1a1 <- glmmTMB(Beetle.Wgt~1, data=dd) # null model
lm1a2 <- glmmTMB(Beetle.Wgt~Treatment, data=dd) # linear model
lm1a3 <- glmmTMB(Beetle.Wgt~Treatment+I(Treatment^2), data=dd) # quadratic model 
anova(lm1a1,lm1a2,lm1a3)

Anova(lm1a2)
summary(lm1a2)

dd1ae <- ggplot(dd%>% filter(Beetle.Wgt<20), aes(x=Treatment, y=Beetle.Wgt))+
  geom_smooth(method="lm", formula = y~x+I(x^2))+
  geom_point() +
  labs(y="Change in beetles mass (mg)", x="Diosgenin in diet (mg/g)")+
  theme_bw(base_size = 24)

lm2a1 <- glmmTMB(Frass.Wgt~1, data=dd)
lm2a2 <- glmmTMB(Frass.Wgt~Treatment, data=dd)
lm2a3 <- glmmTMB(Frass.Wgt~Treatment+I(Treatment^2), data=dd)
anova(lm2a1,lm2a2,lm2a3)

Anova(lm2a2)
MuMIn::r.squaredGLMM(lm2a2)

dd2be <- ggplot(dd , aes(x=Treatment*100, y=Frass.Wgt*1000, fill=Treatment))+
  geom_smooth(method="lm", formula = y~x, color="#395D9C", fill="#395D9C", lwd=2)+
  geom_jitter(height=0, width=.1, pch=21, size=3) +
  scale_fill_viridis(direction=-1, option="G")+
  labs(y="Frass mass (mg)", x="Diosgenin in diet (mg/g)")+
  theme_bw(base_size = 24)+
  theme(legend.position = "none")



ggsave("ddPlot2egena.tiff", dd2be, width=3, height=4, units="in", dpi=600, compression = "lzw")



# bulbil feeding data ####
bd<-read.csv("Copy of Bulbil_FeedingExp_01_04_24.csv")
head(bd)

bm1 <- glmmTMB(Frass_weight~Treatment, data=bd )
Anova(bm1)
emmeans(bm1, pairwise~Treatment)

bd1 <- ggplot(bd , aes(x=Treatment, y=log(Frass_weight*1000)))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(height = 0, width=.2)+
  labs(y="Frass mass (mg)", x="Treatment")+
  theme_bw(base_size = 24)

bm2 <- glmmTMB(Pupa_weight..g.~Treatment, data=bd )
Anova(bm2)
emmeans(bm2, pairwise~Treatment)

bd2 <- ggplot(bd , aes(x=Treatment, y=Pupa_weight..g.*1000))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(height = 0, width=.2)+
  labs(y="Pupae weight (mg)", x="Treatment")+
  theme_bw(base_size = 24)


bm3 <- glmmTMB(Tunnel_weight..g.~Treatment, data=bd )
Anova(bm3)
emmeans(bm3, pairwise~Treatment)

bd3 <- ggplot(bd , aes(x=Treatment, y=Tunnel_weight..g.))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(height = 0, width=.2)+
  labs(y="Weight of wax", x="Treatment")+
  theme_bw(base_size = 24)


bm4 <- glmmTMB(Pupa_weight..g.*1000~Tunnel_weight..g.*Treatment, data=bd)
Anova(bm4)
emmeans(bm4, pairwise~Treatment)
emtrends(bm4, ~Treatment, var="Tunnel_weight..g.", infer=T)

bd4 <- ggplot(bd , aes(x=Tunnel_weight..g., y=Pupa_weight..g.*1000))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~Treatment)+
  labs(y="Pupae weight (mg)", x="Weight of wax")+
  theme_bw(base_size = 24)

bdplot <- bd3+bd2
ggsave("bdPlot.tiff", bdplot, width=15, height=6, units="in", dpi=600, compression = "lzw")
