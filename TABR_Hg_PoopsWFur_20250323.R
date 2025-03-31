#---
#title: TABR_Hg_PoopsWFur
#author: Molly Simonis
#date: 2024-03-23
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
library(geomtextpath)

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

poop_fur_Hg$sample_type[poop_fur_Hg$sample_type == 'feces']<- 'Feces'
poop_fur_Hg$sample_type[poop_fur_Hg$sample_type == 'fur']<- 'Fur'


#Make current df wide for sample_type for descriptive results
poop_fur_Hg_wide<- dcast(poop_fur_Hg, UID + age + date_collected + mass + sex + reproductive.status ~ sample_type, value.var = 'Hg_mgkg')

#order repstats
poop_fur_Hg$reproductive.status<- factor(poop_fur_Hg$reproductive.status, 
                                         levels = c('nonreproductive', 'pregnant',
                                                    'lactating', 'postlactating'))

poop_fur_Hg_wide$reproductive.status<- factor(poop_fur_Hg_wide$reproductive.status, 
                                              levels = c('nonreproductive', 'pregnant',
                                                         'lactating', 'postlactating'))



#descriptive statistics
#average Hg
round(mean(poop_fur_Hg$Hg_mgkg), 2) #0.66

#std error Hg
round(sd(poop_fur_Hg$Hg_mgkg)/sqrt(41), 2) #0.09

#number of bats each date
table(poop_fur_Hg_wide$date_collected)
#18-Jul-23 18-Jun-23 19-Jun-23   unknown 
#.      22        12        13         1 

#number of repstats
table(poop_fur_Hg_wide$reproductive.status[poop_fur_Hg_wide$sex == 'female'])
#nonreproductive        pregnant       lactating   postlactating 
#             5               9              25               5 
table(poop_fur_Hg_wide$reproductive.status[poop_fur_Hg_wide$sex == 'male'])
#nonreproductive 
#     3


#number of ages
table(poop_fur_Hg_wide$age[poop_fur_Hg_wide$sex == 'female'])
#adult juvenile subadult 
#   41        2        1
table(poop_fur_Hg_wide$age[poop_fur_Hg_wide$sex == 'male'])
#juvenile 
#       3

#raw means and ranges
round(mean(poop_fur_Hg$Hg_mgkg[poop_fur_Hg$sample_type == 'Feces']), 
      2) #0.19
round(sd(poop_fur_Hg$Hg_mgkg[poop_fur_Hg$sample_type == 'Feces'])/sqrt(48), 
      2) #0.01
round(min(poop_fur_Hg$Hg_mgkg[poop_fur_Hg$sample_type == 'Feces']), 
      2) #0.10
round(max(poop_fur_Hg$Hg_mgkg[poop_fur_Hg$sample_type == 'Feces']), 
      2) #0.43

round(mean(poop_fur_Hg$Hg_mgkg[poop_fur_Hg$sample_type == 'Fur']), 
      2) #1.14
round(sd(poop_fur_Hg$Hg_mgkg[poop_fur_Hg$sample_type == 'Fur'])/sqrt(48), 
      2) #0.07
round(min(poop_fur_Hg$Hg_mgkg[poop_fur_Hg$sample_type == 'Fur']), 
      2) #0.39
round(max(poop_fur_Hg$Hg_mgkg[poop_fur_Hg$sample_type == 'Fur']), 
      2) #2.12



#create model for differences between pooled fecal and fur Hg
#Using both UID and date_collected as random effects results in singular fit
#Using only UID as random effect instead
Hg_glmer<- glmer(Hg_mgkg ~ sample_type + (1|UID) + (1|date_collected), data = poop_fur_Hg, family = Gamma(link = 'log'))

Anova(Hg_glmer, type = 2)
#Chisq = 764.12, fe df = 1, re df = 47, P < 0.0001

summary(Hg_glmer)
#all looks ok

round(r.squaredGLMM(Hg_glmer), 2)
#R2m = 0.83, R2c = 0.88

#check assumptions
plot_model(Hg_glmer, type = 'diag')
#looks as good as it'll get

#get mean values for sample type only since no sig interx
poop_fur_Hg_means<- data.frame(emmeans(Hg_glmer, ~sample_type, type = "response"))
names(poop_fur_Hg_means)<- c('sample_type', 'Hg_mgkg', 'SE', 'df', 'lowerCI', 'upperCI')

round(poop_fur_Hg_means$Hg_mgkg[poop_fur_Hg_means$sample_type == 'Fur'], 2) /
        round(poop_fur_Hg_means$Hg_mgkg[poop_fur_Hg_means$sample_type == 'Feces'], 2)
#fur is 6.06 times greater than feces

round(poop_fur_Hg_means$Hg_mgkg[poop_fur_Hg_means$sample_type == 'Fur'], 2)
round(poop_fur_Hg_means$lowerCI[poop_fur_Hg_means$sample_type == 'Fur'], 2)
round(poop_fur_Hg_means$upperCI[poop_fur_Hg_means$sample_type == 'Fur'], 2)
#1.09 (0.95, 1.25) for mean and CIs for fur

round(poop_fur_Hg_means$Hg_mgkg[poop_fur_Hg_means$sample_type == 'Feces'], 2)
round(poop_fur_Hg_means$lowerCI[poop_fur_Hg_means$sample_type == 'Feces'], 2)
round(poop_fur_Hg_means$upperCI[poop_fur_Hg_means$sample_type == 'Feces'], 2)
#0.18 (0.16, 0.20) for mean and CIs for feces

#CI's are 5.94 (LCL) and 6.25 (UCL) times greater in fur than feces
#suggested conversion factor = 6 [5.94, 6.25]

confint(contrast(emmeans(Hg_glmer, ~sample_type, type = "response"), method = "pairwise"))
#contrast    ratio     SE  df asymp.LCL asymp.UCL
#feces / fur 0.164  0.0107 Inf    0.144     0.186
#
#Confidence level used: 0.95 
#Intervals are back-transformed from the log scale 



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
  labs(x= "", y= expression(paste("THg Concentration (", mu, "g/g)")))



#individual THg relationships
Hg_glmer2<- glmer(Fur ~ Feces + (1|date_collected), data = poop_fur_Hg_wide, family = Gamma(link = 'log'))

Anova(Hg_glmer2, type = 2)

summary(Hg_glmer2)

r.squaredGLMM(Hg_glmer2)

plot_model(Hg_glmer2, type = 'diag')
#as good as it will get

#determine residuals of data points from a 6:1 line
#use equation resid = y - ax - b
#so resid from a 6:1 line = fur - 6*feces - 0
resid_6_line<- poop_fur_Hg_wide$Fur - 6*poop_fur_Hg_wide$Feces

table(resid_6_line<0)
#28 above the line, 58%
#20 below the line, 42%
#matches to a point count on Fig 2

max(resid_6_line) #largest resid above the line = 1.09
min(resid_6_line) #largest resid below the line = 2.15
min(abs(resid_6_line)) #point closest to line (sits positive in data)
max(abs(resid_6_line)) #point furthest from line (sits negative in data, same as min resid above)

mean(resid_6_line)
round(sd(resid_6_line)/sqrt(48), 2)



#Plot Fig 2
ab_lab<- 'y = 6x'

ab_CI_df<- as.data.frame(cbind(Feces = seq(0, 0.5, by = 0.005), 
                               Fur = 6*seq(0, 0.5, by = 0.005), 
                               LCL = 6*seq(0, 0.5, by = 0.005) - 0.06, 
                               UCL = 6*seq(0, 0.5, by = 0.005) + 0.25))


ggplot(data = poop_fur_Hg_wide, aes(x = Feces, y = Fur)) + 
  geom_point(size = 3, alpha = 0.75, color = '#28788EFF') + 
  geom_line(data = ab_CI_df, aes(x = Feces, y = Fur),
            color = 'gray35', lwd = 1.5, linetype = 'dashed') +
  geom_ribbon(data = ab_CI_df, aes(ymin = Fur - 0.06, ymax = Fur + 0.25),
              fill = 'gray35', alpha = 0.2) +
  coord_cartesian(xlim = c(min(poop_fur_Hg_wide$Feces), max(poop_fur_Hg_wide$Feces)), 
                  ylim = c(min(poop_fur_Hg_wide$Fur), max(poop_fur_Hg_wide$Fur))) +
  annotate(geom = "text", x = 0.25, y = 1.25, size = 5.5,
           label = ab_lab, color = "gray35", angle = 42, vjust = -3.5, hjust = -0.25) +
  th +
  labs(x = expression(paste("Fecal THg Concentration (", mu, "g/g)")), 
       y = expression(paste("Fur THg Concentration (", mu, "g/g)")))

