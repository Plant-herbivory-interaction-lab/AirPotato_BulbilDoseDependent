library(tidyverse)
library(car)
library(emmeans)
library(glmmTMB)
library(patchwork)
library(GGally)

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
db <- left_join(dd,bb) %>% 
  filter(Use=="Feeding", Bulbil_category!='Extra small')

# explore data ####
head(db)

ggplot(db, aes(x=Treatment, y=Feeding_slice)) +
  geom_boxplot()+
  facet_wrap(~Bulbil_category) ## Extra smalls are very small, so remove 



## Tunnel and pupae weight ####
ggplot(db, aes(x=Treatment, y=Tunnel_weight..g.)) +
  geom_boxplot()
ggplot(db, aes(x=Treatment, y=Pupa_weight..g.)) +
  geom_boxplot()
ggplot(db, aes(x=Treatment, y=Frass_weight)) +
  geom_boxplot()
ggplot(db, aes(x=Treatment, y=Emerged_beetle_weight)) +
  geom_boxplot()

db_resp <- db %>% select(Treatment, Tunnel_weight..g., Pupa_weight..g., Frass_weight, Emerged_beetle_weight)
ggpairs(db_resp)

# analyze data ####

mod1 <- glmmTMB(Tunnel_weight..g. ~ Treatment, data=db)
mod1 <- glmmTMB(Emerged_beetle_weight ~ Treatment, data=db)
Anova(mod1)

# Pupa weight ####
ggplot(dd, aes(x=Treatment, y=Pupa_weight..g.)) +
  geom_boxplot() 

mod2a <- glmmTMB((Pupa_weight..g.) ~ Treatment, data=dd)
Anova(mod2a)
DHARMa::simulateResiduals(mod2a, plot=T)
performance::r2(mod2a)
emmeans(mod2a, pairwise ~ Treatment)

mod2 <- glmmTMB(Pupa_weight..g. ~ Treatment * Tunnel_weight..g., data=dd)
Anova(mod2)
emmeans(mod2, pairwise ~ Treatment)
emtrends(mod2, var="Tunnel_weight..g.", pairwise~Treatment, infer=T)

# Frass weight ####
ggplot(dd %>% filter(Frass_weight>0), aes(x=Treatment, y=Frass_weight)) +
  geom_boxplot()

mod3 <- glmmTMB(Frass_weight ~ Treatment, data=dd %>% filter(Frass_weight>0))
Anova(mod3)

ggplot(dd %>% filter(Frass_weight>0), aes(x=Tunnel_weight..g., y=Pupa_weight..g.)) +
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~Treatment)+
  theme_bw(base_size = 18)

ggplot(dd %>% filter(Frass_weight>0), aes(x=Frass_weight, y=Pupa_weight..g.)) +
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~Bulbil_category)+
  theme_bw(base_size = 18)

ggplot(dd %>% filter(Frass_weight>0), aes(x=Frass_weight, y=Pupa_weight..g.)) +
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~Treatment)+
  theme_bw(base_size = 18)

plot1 <- ggplot(dd, aes(x=Treatment, y=Tunnel_weight..g.)) +
  geom_boxplot(outlier.shape=NA, fill="firebrick1", lwd=1.25) +
  geom_jitter(height=0, width=.2) +
  ylab("Amount of bulbil consumed (g)")+
  theme_bw(base_size = 24)

plot2 <- ggplot(dd, aes(x=Treatment, y=Pupa_weight..g.*1000)) +
  geom_boxplot(outlier.shape=NA, fill="firebrick1", lwd=1.25) +
  geom_jitter(height=0, width=.2) +
  ylab("Pupa weight (mg)")+
  theme_bw(base_size = 24)

beetle_plot <- plot1 + plot2
ggsave("beetle_plot.tiff", beetle_plot, width=16, height=6, units="in", dpi=600, compression = "lzw")

