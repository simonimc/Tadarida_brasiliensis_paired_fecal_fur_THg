#---
#title: TABR_Hg_PoopsWFur
#author: Molly Simonis
#date: 2024-03-16
#---

#This R code is for understanding relationships between paired samples of fecal and fur total Hg,
#collected from T. brasiliensis in summer 2023. 

#Make sure everything is clear in the workspace 
rm(list=ls())##clear environment
graphics.off()##clear graphics
gc()##free up some RAM

#turn on packages
library(ggplot2)
library(lmerTest)
library(lme4)
library(car)
library(sjPlot)
library(glmmTMB)
library(MuMIn)
library(emmeans)
library(reshape2)

#make general theme for plots later
th=theme_bw()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  theme(plot.title = element_text(hjust = 0.5))

#read in data
poop_fur_Hg<- read.csv('poop_fur_Hg.csv', sep = ',', header = T)

#Make current df wide for sample_type for later analyses
poop_fur_Hg_wide<- dcast(poop_fur_Hg, UID + date_collected + mass + sex + age + reproductive.status ~ sample_type, value.var = 'Hg_mgkg')


#descriptive statistics
#average Hg
mean(poop_fur_Hg$Hg_mgkg)
#std error Hg
sd(poop_fur_Hg$Hg_mgkg)/sqrt(48)
#number of males and females in sample
table(poop_fur_Hg_wide$sex)
#number of bats each date
table(poop_fur_Hg_wide$date_collected)
#number of repstats
table(poop_fur_Hg_wide$reproductive.status)
table(poop_fur_Hg_wide$reproductive.status[poop_fur_Hg_wide$sex == 'male'])
#number of ages
table(poop_fur_Hg_wide$age)
table(poop_fur_Hg_wide$age[poop_fur_Hg_wide$sex == 'male'])
table(poop_fur_Hg_wide$age[poop_fur_Hg_wide$sex == 'female'])

mean(poop_fur_Hg$Hg_mgkg[poop_fur_Hg$sample_type == 'feces']) #0.1851058
sd(poop_fur_Hg$Hg_mgkg[poop_fur_Hg$sample_type == 'feces'])/sqrt(48) #0.01021918

mean(poop_fur_Hg$Hg_mgkg[poop_fur_Hg$sample_type == 'fur']) #1.142281
sd(poop_fur_Hg$Hg_mgkg[poop_fur_Hg$sample_type == 'fur'])/sqrt(48) #0.06816821

#create model for differences between fecal and fur Hg
Hg_glmer<- glmer(Hg_mgkg ~ sample_type + (1|UID) + (1|date_collected), data = poop_fur_Hg, family = Gamma(link = 'log'))

Anova(Hg_glmer, type = 2)
#Analysis of Deviance Table (Type II Wald chisquare tests)
#
#Response: Hg_mgkg
#Chisq Df Pr(>Chisq)    
#sample_type 764.12  1  < 2.2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(Hg_glmer)
#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#Family: Gamma  ( log )
#Formula: Hg_mgkg ~ sample_type + (1 | UID) + (1 | date_collected)
#Data: poop_fur_Hg
#
#AIC      BIC   logLik deviance df.resid 
#-68.5    -55.7     39.3    -78.5       91 
#
#Scaled residuals: 
#  Min       1Q   Median       3Q      Max 
#-1.92493 -0.71716 -0.08268  0.71236  3.07681 
#
#Random effects:
#  Groups         Name        Variance  Std.Dev.
#UID            (Intercept) 0.0493437 0.2221  
#date_collected (Intercept) 0.0005903 0.0243  
#Residual                   0.1182630 0.3439  
#Number of obs: 96, groups:  UID, 48; date_collected, 4
#
#Fixed effects:
#  Estimate Std. Error t value Pr(>|z|)    
#(Intercept)    -1.72557    0.07067  -24.42   <2e-16 ***
#  sample_typefur  1.80921    0.06545   27.64   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Correlation of Fixed Effects:
#  (Intr)
#sampl_typfr -0.471

r.squaredGLMM(Hg_glmer)
#           R2m       R2c
#delta     0.8309789 0.8811575
#lognormal 0.8364309 0.8869388
#trigamma  0.8249537 0.8747686

#check assumptions
plot_model(Hg_glmer, type = 'diag')

#get mean values
poop_fur_Hg_means<- data.frame(emmeans(Hg_glmer, ~sample_type, type = "response"))
names(poop_fur_Hg_means)<- c('sample_type', 'Hg_mgkg', 'SE', 'df', 'lowerCI', 'upperCI')

1.0872378/0.1780711 #fur is 6.11 greater than feces
#feces  (0.18 [0.16, 0.21])  
#fur    (1.087 [0.95, 1.25])


confint(contrast(emmeans(Hg_glmer, ~sample_type, type = "response"), method = "pairwise"))
#contrast    ratio     SE  df asymp.LCL asymp.UCL
#feces / fur 0.164 0.0107 Inf     0.144     0.186
#
#Confidence level used: 0.95 
#Intervals are back-transformed from the log scale 

#backtransform estimates for mean and CI's of the contrast
log(0.164) #estimate = -1.807889
log(0.144) #LCL = -1.937942
log(0.186) #UCL = -1.682009
#the mean difference between fecal and fur THg is that feces is 1.81 mg/kg
#less than fur on average (-1.81 [-1.94, -1.68])
#
#plot Fig 1
ggplot(data = poop_fur_Hg, aes(x = sample_type, y = Hg_mgkg)) +
  geom_path(aes(group = UID, color = UID), alpha = 0.30, lwd = 0.75) +
  geom_point(aes(group = UID, color = UID), alpha = 0.30) +
  scale_color_viridis_d()+
  geom_path(data = poop_fur_Hg_means, aes(x = sample_type, y = Hg_mgkg, group = 1), color = 'black', linewidth = 2)+
  geom_pointrange(data = poop_fur_Hg_means, aes(x = sample_type, y = Hg_mgkg, ymin = lowerCI, ymax = upperCI), 
                  size = 1, lwd = 1, shape = 22, fill = "black") +
  ylim(0, max(poop_fur_Hg$Hg_mgkg))+
  th + theme(legend.position = "none") + #theme(plot.margin = unit(c(0, 0, 0, 0), 'mm')) +
  labs(x="Sample Matrix", y="THg Concentration (mg/kg)")


#create model for a direct relationship between fecal and fur Hg
#create model
Hg_glmer2<- glmer(fur ~ feces + (1|date_collected), data = poop_fur_Hg_wide, family = Gamma(link = 'log'))
#transforming data doesn't really give anything better or worse

Anova(Hg_glmer2, type = 2)
#Analysis of Deviance Table (Type II Wald chisquare tests)
#
#Response: fur
#Chisq Df Pr(>Chisq)
#feces 0.2577  1     0.6117

summary(Hg_glmer2)
#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#Family: Gamma  ( log )
#Formula: fur ~ feces + (1 | date_collected)
#Data: poop_fur_Hg_wide
#
#AIC      BIC   logLik deviance df.resid 
#62.5     70.0    -27.3     54.5       44 
#
#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-1.6193 -0.9266 -0.1938  0.8872  1.7100 
#
#Random effects:
#  Groups         Name        Variance Std.Dev.
#date_collected (Intercept) 0.02155  0.1468  
#Residual                   0.13694  0.3701  
#Number of obs: 48, groups:  date_collected, 4
#
#Fixed effects:
#  Estimate Std. Error t value Pr(>|z|)
#(Intercept)   0.1502     0.2186   0.687    0.492
#feces         0.4412     0.8691   0.508    0.612
#
#Correlation of Fixed Effects:
#  (Intr)
#feces -0.677

r.squaredGLMM(Hg_glmer2)
#            R2m       R2c
#delta     0.006119469 0.1412590
#lognormal 0.006468289 0.1493110
#trigamma  0.005765076 0.1330783


#plot Fig 2
ggplot(data = poop_fur_Hg_wide, aes(x = feces, y = fur)) +
  geom_point(color = 'turquoise3') +
  th + 
  labs(x="Fecal THg Concentration (mg/kg)", y="Fur THg Concentration (mg/kg)") +
  ylim(0, max(poop_fur_Hg_wide$fur))+
  xlim(0, max(poop_fur_Hg_wide$feces))




#Hg with mass
Hg_glmer3<- glmer(mass ~ Hg_mgkg*sample_type + (1|UID) + (1|date_collected), data = poop_fur_Hg, family = Gamma(link = 'log'))
#using site as random effect screws up model convergence

Anova(Hg_glmer3, type = 2)
#Analysis of Deviance Table (Type II Wald chisquare tests)
#
#Response: mass
#Chisq Df Pr(>Chisq)
#Hg_mgkg                 0  1     1.0000
#sample_type             0  1     1.0000
#Hg_mgkg:sample_type     0  1     0.9999

summary(Hg_glmer3)
#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#Family: Gamma  ( log )
#Formula: mass ~ Hg_mgkg * sample_type + (1 | UID) + (1 | date_collected)
#Data: poop_fur_Hg
#
#AIC      BIC   logLik deviance df.resid 
#-2002.4  -1984.9   1008.2  -2016.4       83 
#
#Scaled residuals: 
#  Min         1Q     Median         3Q        Max 
#-3.727e-04 -9.130e-05  1.022e-05  8.920e-05  3.900e-04 
#
#Random effects:
#  Groups         Name        Variance  Std.Dev. 
#UID            (Intercept) 1.055e-02 1.027e-01
#date_collected (Intercept) 1.149e-14 1.072e-07
#Residual                   5.008e-10 2.238e-05
#Number of obs: 90, groups:  UID, 45; date_collected, 3
#
#Fixed effects:
#  Estimate Std. Error t value Pr(>|z|)    
#(Intercept)             2.467e+00  1.531e-02   161.2   <2e-16 ***
#  Hg_mgkg                 9.124e-09  7.038e-05     0.0        1    
#sample_typefur          1.303e-09  1.796e-05     0.0        1    
#Hg_mgkg:sample_typefur -8.962e-09  7.123e-05     0.0        1    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Correlation of Fixed Effects:
#  (Intr) Hg_mgk smpl_t
#Hg_mgkg     -0.001              
#sampl_typfr -0.001  0.717       
#Hg_mgkg:sm_  0.001 -0.989 -0.804
#optimizer (Nelder_Mead) convergence code: 0 (OK)
#Gradient contains NAs

r.squaredGLMM(Hg_glmer3)
#                   R2m R2c
#delta     1.894612e-17   1
#lognormal 1.894612e-17   1
#trigamma  1.894612e-17   1

#plot Fig 3
ggplot(data = poop_fur_Hg, aes(x = mass, y = Hg_mgkg, color = sample_type)) +
  geom_point() +
  scale_color_viridis_d(direction = -1, begin = 0, end =  0.55)+
  ylim(0, max(poop_fur_Hg$Hg_mgkg))+
  th + 
  theme(legend.position = 'top') + theme(legend.text = element_text(size = 12)) +
  labs(x=expression(paste(italic('T. brasiliensis'), ' Body Mass (g)')), 
       y="THg Concentration (mg/kg)", color = 'Sample Matrix', fill = 'Sample Matrix') 
