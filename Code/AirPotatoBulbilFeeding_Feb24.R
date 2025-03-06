library(tidyverse)
library(car)
library(emmeans)
library(glmmTMB)
library(patchwork)
library(GGally)
library(skimr)
library(DHARMa)

# LOAD DATA FILES ####

## load and process bulbil data file ####
bb<-read.csv("Data/Bulbil_FeedingExp_04_08_24_bulbil.csv") #%>% #filter(Wave!='1') #%>%
  # mutate(Treatment = fct_recode(as.factor(Treatment), "L. cheni herbivory" = "Beetle ")) %>% 
  # mutate(Treatment = fct_relevel(as.factor(Treatment), "Control ","L. cheni herbivory","JA", "SA"))
head(bb)  

## load and process feeding file ####
dd<-read.csv("Data/Bulbil_FeedingExp_04_08_24_feeding.csv") %>% #filter(Wave!='1') %>%
  # mutate(Treatment = fct_recode(as.factor(Treatment), "L. cheni herbivory" = "Beetle ")) %>% 
  # mutate(Treatment = fct_relevel(as.factor(Treatment), "Control ","L. cheni herbivory","JA", "SA")) %>% 
  drop_na(Tunnel_weight..g.) 
head(dd)  

### join feeding and bulbil files ####
## also remove spaces in Beetle and Control and reorder levels
db <- left_join(dd,bb) %>% 
  filter(Use=="Feeding", Bulbil_category!='Extra small') %>% 
  mutate(Treatment=gsub("Beetle ", "Beetle",Treatment),
         Treatment=gsub("Control ", "Control",Treatment),
         Treatment = factor(Treatment, levels=c("Control", "Beetle", "JA", "SA")))
  

# explore data ####
head(db)

ggplot(db, aes(x=Treatment, y=Feeding_slice)) +
  geom_boxplot()+
  facet_wrap(~Bulbil_category) ## Extra smalls are very small, so remove 

ggplot(db, aes(x=Treatment, y=Frass_weight)) +
  geom_boxplot()+
  facet_wrap(~Wave)

## Tunnel and pupae weight ####
ggplot(db, aes(x=Treatment, y=Tunnel_weight..g.)) +
  geom_boxplot()
ggplot(db, aes(x=Treatment, y=Pupa_weight..g.)) +
  geom_boxplot()
ggplot(db, aes(x=Treatment, y=Frass_weight)) +
  geom_boxplot()
ggplot(db, aes(x=Treatment, y=Emerged_beetle_weight)) +
  geom_boxplot()

db_resp <- db %>% select(Treatment, Tunnel_weight..g., Pupa_weight..g., Frass_weight, Emerged_beetle_weight, Feeding_slice)
summary(db_resp)
skim(db_resp)
ggpairs(db_resp)

## REPORT THESE NOTES IN METHODS
## Bulbil categories doesn't seem to matter much, once Extra Smalls were excluded
## Size of Feeding_slice doesn't seem to matter much
## Not many moralities (n = 13 of 99 trials) and they are spread across treatments, probably ok to ignore them

# analyze data ####
mm1 <- lm(cbind(Tunnel_weight..g., Frass_weight, Pupa_weight..g., Emerged_beetle_weight) ~ Treatment, data=db)
Manova(mm1, test.statistic="Wilks")  ## report results for manova, approx. F and p-value

## indivdual Anova's ####
m1 <- glmmTMB(Tunnel_weight..g. ~ Treatment , data=db)
Anova(m1)
simulateResiduals(m1, plot=T)

m2 <- glmmTMB(Frass_weight ~ Treatment #+ Feeding_slice
              , data=db)
Anova(m2)
simulateResiduals(m2, plot=T)

m3 <- glmmTMB(Pupa_weight..g. ~ Treatment, data=db)
Anova(m3)
emmeans(m3, pairwise ~ Treatment)
simulateResiduals(m3, plot=T)

m4 <- glmmTMB(Emerged_beetle_weight ~ Treatment, data=db)
Anova(m4)
emmeans(m4, pairwise ~ Treatment)
simulateResiduals(m4, plot=T)


# make figures ####
bd1 <- ggplot(db , aes(x=Treatment, y=Tunnel_weight..g.*1000))+  
  geom_boxplot(alpha = 0.4, outlier.shape = NA, aes(fill = Treatment))+  
  geom_jitter(height = 0, width = 0.1, size= 2, alpha=0.7)+  
  labs(y="Wax weight (mg)", x="")+
  theme_bw(base_size = 24)

bd2 <- ggplot(db , aes(x=Treatment, y=Frass_weight*1000))+  
  geom_boxplot(alpha = 0.4, outlier.shape = NA, aes(fill = Treatment))+  
  geom_jitter(height = 0, width = 0.1, size= 2, alpha=0.7)+  
  labs(y="Frass weight (mg)", x="")+
  theme_bw(base_size = 24)

bd3 <- ggplot(db , aes(x=Treatment, y=Pupa_weight..g.*1000))+  
  geom_boxplot(alpha = 0.4, outlier.shape = NA, aes(fill = Treatment))+  
  geom_jitter(height = 0, width = 0.1, size= 2, alpha=0.7)+  
  labs(y="Pupae weight (mg)", x="")+
  theme_bw(base_size = 24)

bd4 <- ggplot(db , aes(x=Treatment, y=Emerged_beetle_weight*1000))+  
  geom_boxplot(alpha = 0.4, outlier.shape = NA, aes(fill = Treatment))+  
  geom_jitter(height = 0, width = 0.1, size= 2, alpha=0.7)+  
  labs(y="Beetle weight (mg)", x="")+
  theme_bw(base_size = 24)

beetle_plot <- (bd1+bd2)/(bd3+bd4) +
  plot_annotation(tag_levels = 'a') + plot_layout(guides='collect') & theme(legend.position='none')

ggsave("beetle_plot.tiff", beetle_plot, width=12, height=9, units="in", dpi=600, compression = "lzw", path="Outputs")

## Figures for SEB March 2025 presentation ####
bd2s <- ggplot(db , aes(x=Treatment, y=Frass_weight*1000))+  
  geom_boxplot(alpha = 0.4, outlier.shape = NA, fill = "Red")+  
  geom_jitter(height = 0, width = 0.1, size= 3, alpha=0.7)+  
  labs(y="Frass weight (mg)", x="Treatment")+
  theme_bw(base_size = 24)
bd3s <- ggplot(db , aes(x=Treatment, y=Pupa_weight..g.*1000))+  
  geom_boxplot(alpha = 0.4, outlier.shape = NA, fill = "Red")+  
  geom_jitter(height = 0, width = 0.1, size= 3, alpha=0.7)+  
  labs(y="Pupae weight (mg)", x="Treatment")+
  theme_bw(base_size = 24)

beetle_plotSEB <- (bd2s+bd3s)+
  plot_annotation(tag_levels = 'a') + plot_layout(guides='collect') & theme(legend.position='none')

ggsave("beetle_plotSEB.tiff", beetle_plotSEB, width=15, height=6, units="in", dpi=600, compression = "lzw", path="Outputs")



# EXTRA analysis for growth rate -- correlations between wax or frass weight and pupae weight ####
## other figures ####
mod1 <- glmmTMB(Pupa_weight..g.*1000~Tunnel_weight..g.*Treatment, data=db)
Anova(mod1)
emmeans(mod1, pairwise~Treatment)
emtrends(mod1, ~Treatment, var="Tunnel_weight..g.", infer=T)

mod2 <- glmmTMB(Pupa_weight..g.*1000~Frass_weight*Treatment, data=db)
Anova(mod2)
emmeans(mod2, pairwise~Treatment)
emtrends(mod2, ~Treatment, var="Frass_weight", infer=T)

# can color these graphs the same as the boxplot if wanted
extraplot1 <- ggplot(db , aes(x=Tunnel_weight..g., y=Pupa_weight..g.*1000))+
  geom_point()+
  geom_smooth(method="lm", color="black", fill="red")+
  facet_wrap(~Treatment)+
  labs(y="Pupae weight (mg)", x="Weight of wax (g)")+
  theme_bw(base_size = 24)

extraplot2 <- ggplot(db , aes(x=Frass_weight, y=Pupa_weight..g.*1000))+
  geom_point(size=2)+
  geom_smooth(method="lm", color="black", fill="red")+
  facet_wrap(~Treatment)+
  labs(y="Pupae weight (mg)", x="Frass weight (g)")+
  theme_bw(base_size = 24)

extraplot1 <- ggplot(db, aes(x = Tunnel_weight..g., y = Pupa_weight..g.*1000, color = Treatment)) +
  geom_point(alpha = 0.7, size = 2) + 
  geom_smooth(method = "lm", color = "black", aes(fill = Treatment), alpha = 0.2) + 
  facet_wrap(~Treatment, labeller = labeller(Treatment = c(
    "Control" = "a.  Control", 
    "Beetle" = "b.  Beetle", 
    "JA" = "c.  Jasmonic Acid", 
    "SA" = "d.  Salicylic Acid"))) + 
  labs(y = "Pupae weight (mg)", x = "Weight of wax (g)") +
  theme_bw(base_size = 24) +
  theme( strip.text = element_text(hjust = 0, face = "bold", size = 24),
         legend.position = "none");extraplot1


extraplot2 <- ggplot(db, aes(x = Frass_weight, y = Pupa_weight..g.*1000, color = Treatment)) +
  geom_point(alpha = 0.7, size = 2) + 
  geom_smooth(method = "lm", color = "black", aes(fill = Treatment), alpha = 0.2) + 
  facet_wrap(~Treatment, labeller = labeller(Treatment = c(
    "Control" = "a.  Control", 
    "Beetle" = "b.  Beetle", 
    "JA" = "c.  Jasmonic Acid", 
    "SA" = "d.  Salicylic Acid"))) + 
  labs(y = "Pupae weight (mg)", x = "Frass weight (mg)") +
  theme_bw(base_size = 24) +
  theme( strip.text = element_text(hjust = 0, face = "bold", size = 24),
         legend.position = "none"); extraplot2



ggsave("extraplot2.tiff", extraplot2, width=12, height=9, units="in", dpi=600, compression = "lzw", path="Outputs")

ggsave("beetle_plot2SEB.tiff", extraplot2, width=10, height=7, units="in", dpi=600, compression = "lzw", path="Outputs")
