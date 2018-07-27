## CSS bird exclusion analysis
# CN
# October 22, 2017
# what is the relationship between plant traits, herbivore communities, and top-down effect of predators?

setwd("/Users/colleennell/Dropbox/Projects/CSS")

library(devtools)
library(ggplot2)
library(dplyr)
library(reshape2)
library(cowplot)
library(tidyr)
library(vegan)

source('/Users/colleennell/Dropbox/rstats/themes/theme_nell.R')

se<-function(x) sd(x, na.rm=TRUE)/sqrt(length(x))

##plant level data on total abundance, plant size, wt, sampling effort
methods<-read.csv('final/CSS_methods_2018.csv')%>%dplyr::select(Sample:dens_g_min)%>%
  mutate(plant_kg = plant_g/1000, plant_vol = pi*(diameter/2)^2*height)
sum(methods$abun_comp) ##2852

#arthropod matrix
arths<-read.csv('final/CSS_master_MAT.csv')
arths[is.na(arths)]<-0
str(arths)
##arthropod groupings to determine pred/herbs
orders<-read.csv('final/CSS_feeding.csv')

arth.mat<-arths%>%dplyr::select(ACAR:WASP)##community matrix
sum(rowSums(arth.mat, na.rm=TRUE)) #2789? which plant 1 or 2) are missing?
##2789 arthropods

##use sizes to calculate approximate arthopod biomass
sizes<-read.csv('final/CSS_arth_sizes.csv')%>%dplyr::select(-L_S, -orig.x)%>%
  transform(Order=as.character(Order), Sample=as.character(Sample))
arth_b<-read.csv('/Users/colleennell/Dropbox/Projects/CAWR/CAWR_arth/data/cawr_arth_biomass.csv')%>%dplyr::select(Order,a,b,bar)
#str(arth_b)
levels(arth_b$Order)<-c('ACAR','ARAN','ARCH','AUCH','CHIL','COLE','DERM','DIPT','ENTO','HEMI','HETE','HYME','ISOP','LEPI',
                        'MANT','NEUR','OPIL','PHAL','ORTH','PSEU','PSOC','RAPH','SIPH','STER','THYS','THYSA')
arth_b$Order<-as.character(arth_b$Order)
sizes$bm_name<-ifelse(sizes$Order == 'COCC', 'COLE', ifelse(sizes$Order == 'MIRID','HETE', ifelse(sizes$Order =='ARCH', 'THYS', sizes$Order)))

##calculate biomass
#Hodar Equation: W=aBL^b
sized<-sizes%>%left_join(arth_b, by=c('bm_name'='Order'))%>%
  mutate(S_bm = (a*1^b)*S, L_bm =(a*4^b)*L, total_bm=S_bm+L_bm)%>%dplyr::select(-bar)%>%transform(treat=as.character(treat))
write.csv(sized, 'final/CSS_arthropod_biomass.csv', row.names=FALSE)
str(sized)
##total arthropo biomass on each plant
arth_bm<-sized%>%group_by(species, Sample, treat)%>%
  summarize(sum=sum(total_bm, na.rm=TRUE))%>%transform(treat=as.character(treat))
##herbivofe biomass
herb_bm<-sized%>%group_by(species, Sample, treat, feed)%>%
  summarize(sum=sum(total_bm, na.rm=TRUE))%>%transform(treat=as.character(treat))%>%
  dcast(species+Sample+treat~feed)
herb_bm[is.na(herb_bm)]<-0
herb_bm$herb_om<-herb_bm$herb+herb_bm$omni
View(herb_bm)
write.csv(herb_bm, 'final/CSS_herbivore_biomass.csv',row.names=FALSE)

##this is a mess 
##add to methods df


######################################
##caluclate biomass per plant
library(lubridate)
str(meth)
##add dry biomass of plant vacuumed in 1 min for each species
##dry biomass vacuumed in 1 min for each plant - arth mg by plant kg
meth<-read.csv('final/CSS_plant_level_data_in.csv')%>%
  dplyr::select(Sample:vaccum, total_abun, plant_g, plant_kg, plant_vol, arth_mg, herb, pred, herb_om)%>%
  mutate(vac_min = as.numeric(seconds(ms(vaccum)))/60, kg_min = plant_kg/vac_min, mean_kg_min=ave(kg_min, sp, FUN=mean))%>%
  mutate(kg_vac = mean_kg_min*vac_min, abun_kg = total_abun/kg_vac, mg_kg = arth_mg/((kg_vac+plant_kg)/2),
         herb_kg=herb/((kg_vac+plant_kg)/2), ho_kg=herb_om/((kg_vac+plant_kg)/2), pred_kg=pred/((kg_vac+plant_kg)/2), 
         po_kg=(pred+omni)/((kg_vac+plant_kg)/2))
str(meth)
write.csv(meth, 'final/CSS_plant_level_data_now.csv',row.names=FALSE)


#arthropod densities
ggplot(meth, aes(interaction(Sample, treat), abun_kg, color=treat))+geom_text(aes(label=Sample))+facet_wrap(~sp, scales='free')
#arthropod biomass
ggplot(meth, aes(interaction(Sample, treat), mg_kg, color=treat))+geom_text(aes(label=Sample))+facet_wrap(~sp, scales='free')

dens.spt<-aov(abun_kg~sp*treat, data=meth)
summary(dens.spt) ##sp, sp*teat
shapiro.test(resid(dens.spt))

#average density bu species * treatmetn
dens.sum<-meth%>%
  group_by(sp, treat)%>%
  summarize(dens_mean = mean(abun_kg, na.rm=TRUE), dens_se = se(abun_kg),
            mg_mean = mean(mg_kg, na.rm=TRUE), mg_se=se(mg_kg), mg_sd=sd(mg_kg, na.rm=TRUE), n=length(unique(Sample)))

##ok gigure- needs different format
ggplot(dens.sum, aes(sp, mg_mean))+
  geom_errorbar(aes(ymin=mg_mean-mg_se, ymax=mg_mean+mg_se, group=treat),width=.15, alpha=.75, size=.8, position=position_dodge(.9))+
  geom_point(aes(shape=treat, fill=treat), size=2.8, stroke=1, position=position_dodge(.9))+
  scale_shape_manual(values=c(24,25))+scale_fill_manual(values=c('black','white'))+
  theme_nell()+labs(x='Species', y='Arthropod density (mg kg-1)')

ggplot(dens.sum, aes(sp, mg_mean))+
  geom_bar(stat='identity',aes(fill=treat), size=.6, alpha=.7,color='black', position=position_dodge(.9))+
  geom_errorbar(aes(ymin=mg_mean-mg_se, ymax=mg_mean+mg_se, group=treat),color='black',width=.15, size=.8, position=position_dodge(.9))+
  scale_shape_manual(values=c(24,25))+scale_fill_manual(values=c('black','white'))+
  theme_nell()+labs(x='Species', y='Arthropod density (mg kg-1)')




#############
# bird effects
library(metafor)
##arthropod density in exclusion (direct defense, random sample of exclusion plants)
set.seed(56)

meth<-read.csv('final/CSS_plant_level_data_now.csv')
meth

##now calculate density by assigned calculation
feed.treat<-meth%>%
  group_by(sp, treat)%>%
  summarize(mean=mean(mg_kg, na.rm=TRUE), sd=sd(mg_kg, na.rm=TRUE), se=se(mg_kg), n=length(unique(Sample)))%>%
  melt(id.vars=c('sp','treat'))%>%
  dcast(sp~treat+variable)
feed.treat

birdalt<-summary(escalc('ROM', m1i=C_mean, m2i=T_mean, sd1i=C_sd, sd2i=T_sd, n1i=C_n, n2i=T_n, data=feed.treat, var.names=c('ID','ID_var'), append=TRUE))

alt.lm<-lm(T_mean~ID, data=birdalt)
summary(aov(alt.lm))
shapiro.test(resid(alt.lm))


#Bird effects + SE (using only half of exclusion plants in claculation of direct defense)
ggplot(birdalt, aes(reorder(sp, ID),ID))+
  geom_errorbar(aes(ymin=ID-sei, ymax=ID+sei), width=.15, alpha=.75, size=.85)+
  geom_point(shape=25, fill='white',stroke=1, size=2.8)+
  geom_hline(yintercept=0, lty='dashed', color='grey')+
  theme_nell()+labs(x='Species', y='Bird Effect\n(LRR control/exclusion)')

##randomly assign groups to calculate LRR
excl<-meth%>%filter(treat == 'T')%>%group_by(sp)
exclu.dd<-sample_n(excl, 4)%>%mutate(group = 'DD')%>%ungroup()#sample n rows
sampn<-unique(exclu.dd$Sample)
exlc.id<-excl%>%filter(!(Sample%in%sampn))%>%mutate(group='ID')%>%ungroup()
cont<-meth%>%filter(treat=='C')%>%mutate(group='control')
exlcusion<-rbind(exclu.dd, exlc.id)
alldata<-rbind(exlcusion, cont)
#View(alldata)
#write.csv(alldata, 'CSS_plants_IDgroups.csv', row.names=FALSE)

lrr.df<-alldata%>%group_by(sp, group)%>%
  summarize(mean=mean(mg_kg, na.rm=TRUE), sd=sd(mg_kg, na.rm=TRUE), se=se(mg_kg), n=length(unique(Sample)))%>%
  melt(id.vars=c('sp','group'))%>%
  dcast(sp~group+variable)
lrr.df

#LRR<-lrr.df%>%left_join(feed.treat, by='sp')
#LRR
#birdz<-summary(escalc('ROM', m1i=control_mean, m2i=DD_mean, sd1i=control_sd, sd2i=DD_sd, n1i=control_n, n2i=DD_n, data=LRR, var.names=c('ID','ID_var'), append=TRUE))
#id.lm<-lm(log(T_mean)~ID, data=birdz)
#summary(aov(id.lm))

#birdfx<-summary(escalc('ROM', m1i=control_mean, m2i=DD_mean, sd1i=control_sd, sd2i=DD_sd, n1i=control_n, n2i=DD_n, data=lrr.df, var.names=c('ID','ID_var'), append=TRUE))

#write.csv(birdfx, 'final/CSS_LRR', row.names=FALSE)
LRR<-read.csv('final/CSS_LRR.csv')
View(LRR)
birdfx<-summary(escalc('ROM', m1i=control_mean, m2i=ID_mean, sd1i=control_sd, sd2i=ID_sd, n1i=control_n, n2i=ID_n, data=LRR, var.names=c('ID','ID_var'), append=TRUE))
birdfx

##Did birds consume more frmo plants with higher arthropod densities?
id.lm<-lm(log(1/DD_mean)~ID, data=birdfx)
summary(aov(id.lm))# P = 0.069
summary(id.lm) ##R2 = .3161
shapiro.test(resid(id.lm))

##trade off in direct and indrect defense? - same as earlier LRR?
ggplot(birdfx, aes(DD_mean,ID))+
  geom_smooth(method='lm', se=F,color='grey', lty='dashed', level=.95)+
  geom_point()+
  geom_hline(yintercept=0, lty='dashed', color='grey')+
  geom_text(aes(label=sp), size=4, hjust=.5, vjust=-1)+
  theme_nell()+labs(x='Plant resistance\n(arthropod density in exclusion)', y='Indirect defense\n(LRR control/exclusion)')+
  scale_x_reverse()+scale_y_reverse()

ggplot(LRR, aes(reorder(sp, ID),ID))+
  geom_errorbar(aes(ymin=ID-sei, ymax=ID+sei), width=.15, alpha=.75, size=.85)+
  geom_point(shape=25, fill='white',stroke=1, size=2.8)+
  geom_hline(yintercept=0, lty='dashed', color='grey')+
  theme_nell()+labs(x='Species', y='Indirect defense\n(LRR control/exclusion)')


birdfx  

write.csv(birdfx,'CSS_newtry.csv',row.names=FALSE)

##with errobars
ggplot(birdfx, aes(DD_mean,ID))+
  geom_errorbar(aes(ymin=ID-sei, ymax=ID+sei), width=0, alpha=.75, size=.75)+
  geom_errorbarh(aes(xmin=DD_mean-DD_se*.7, xmax=DD_mean+DD_se*.7), width=0, alpha=.75, size=.75)+
  geom_smooth(method='lm', se=F, color='darkgrey', lty='dashed')+
  geom_hline(yintercept=0, lty='dashed', color='grey')+
  geom_point(aes(), shape=25, fill='black',stroke=1, size=2.8)+
  theme_nell()+labs(x='Herbivore resistance\n(density in exclusion)', y='Bird effect\n(LRR control/exclusion)')

birdfx
#Bird effects + SE (using only half of exclusion plants in claculation of direct defense)
ggplot(birdfx, aes(reorder(sp, ID),ID))+
  geom_errorbar(aes(ymin=ID-sei, ymax=ID+sei), width=.15, alpha=.75, size=.85)+
  geom_point(shape=25, fill='white',stroke=1, size=2.8)+
  geom_hline(yintercept=0, lty='dashed', color='grey')+
  theme_nell()+labs(x='Species', y='Bird Effect\n(LRR control/exclusion)')

## bird effect should be the difference in biomass?

##HPQ
##does herbivore density relate to HPQ?
hpq<-read.csv("final/CSS_HPQ_SP.csv")
str(hpq)
hpqs<-hpq%>%group_by(plant)%>%
  summarize(hpq=mean(wt, na.rm=TRUE), hpq_se=se(wt))
hpqs

bird.hpq<-birdfx%>%left_join(hpqs, by=c('sp'='plant'))%>%left_join(feed.treat)
bird.hpq

hpq.lmardo<-lm(DD_mean~log(hpq), data=bird.hpq%>%filter(sp!='ARDO'))
hpq.lm<-lm(DD_mean~log(hpq), data=bird.hpq) #P=0.029, R2 = 0.51
summary(aov(hpq.lm))
summary(hpq.lm) #R2 = 0.59, P=0.0155
shapiro.test(resid(hpq.lm))
##non sig if ARDO is not included

bird.hpq$hpq<-ifelse(bird.hpq$sp=='ARDO', 2.6655, bird.hpq$hpq)
bird.hpq$hpq<-ifelse(bird.hpq$sp=='ARCA', 0.21556, bird.hpq$hpq)
bird.hpq$hpq<-ifelse(bird.hpq$sp=='SAME', 0.08756, bird.hpq$hpq)

ggplot(bird.hpq, aes(hpq, DD_mean))+
  geom_smooth(method='lm', se=F, color='grey', lty='dashed')+
  geom_text(aes(label=sp), size=4)+
  theme_nell()+labs(x='Direct defense', y='Host plant quality')+
  scale_x_log10()

##HPQ by the density of herbivores ONLY - things feeding ON THE PLANT
##BIRD EFFECTS - ANY ARTHROPODS (eat all)
##but DD - is just herbivores
#HPQ vs herb density in exclusion
meth<-read.csv('final/CSS_plant_level_data_now.csv')
View(meth)

dd.sp<-meth%>%filter(treat=='T')%>%
  group_by(sp)%>%
  summarize(herb_mean=mean(herb_kg, na.rm=TRUE), herb_se=se(herb_kg),
            ho_mean=mean(ho_kg, na.rm=TRUE), ho_se=se(ho_kg))

hpqq<-bird.hpq%>%left_join(dd.sp, by='sp')
#View(hpqq)
hpqq
hpq.lm<-lm(log(hpq)~ho_mean, data=hpqq)
summary(aov(hpq.lm))#P=0.021
summary(hpq.lm) #R2=.55
shapiro.test(resid(hpq.lm))
str(hpqq)

###good
##HPQ is not correlated with herbivore biomass ine exclusion (P=0.09)
ggplot(hpqq, aes(hpq, ho_mean))+
  geom_smooth(method='lm', se=F, color='grey', lty='dashed')+
  geom_text(aes(label=sp), size=4, position='jitter')+
  theme_nell()+labs(x='Host plant quality', y='Herbivore density')+
  scale_x_log10()

ggplot(hpqq, aes(hpq, ho_mean))+
  geom_errorbar(aes(ymin=ho_mean-ho_se, ymax=ho_mean+ho_se), width=0, alpha=.75, size=.75)+
  geom_errorbarh(aes(xmin=hpq-hpq_se, xmax=hpq+hpq_se), width=0, alpha=.75, size=.75)+
  geom_smooth(method='lm', se=F, color='grey', lty='dashed')+
  geom_jitter(size=2.8, stroke=1, fill='black', shape=24)+
  theme_nell()+labs(x='Host plant quality', y='Herbivore density')+
  scale_x_log10()
hpqq
str(hpqq)
ggplot(hpqq, aes(reorder(sp, hpq), hpq))+geom_bar(stat='identity', fill='grey', color='black')+
  geom_errorbar(aes(ymin=hpq-hpq_se, ymax=hpq+hpq_se), width=.1)+
  theme_nell()+labs(x='Species',y='Larval weight gain (mg)')

summary(aov(hpq~sp, data=hpqq))
###HPQ and ID?
hpq.id<-lm(ID~log(hpq), data=hpqq)
summary(aov(hpq.id))#no
summary(hpq.id) 

ggplot(hpqq, aes(hpq, ID))+
  geom_errorbarh(aes(xmin=hpq-hpq_se, xmax=hpq+hpq_se), width=.1)+
  geom_errorbar(aes(ymin=ID-sei, ymax=ID+sei), width=.1)+
  geom_point(shape=25, size=4, fill='white', stroke=1)+
  theme_nell()+labs(x='HPQ', y='Bird effects')+scale_x_log10()

##differences in predator density with bird exclusion?
preds<-meth%>%
  group_by(sp, treat)%>%
  summarize(pred=mean(pred_kg, na.rm=TRUE), pred_se=se(pred_kg), pro=mean(po_kg, na.rm=TRUE), pro_se=se(po_kg))

comp<-meth%>%group_by(sp)%>%summarize(comp=mean(complexity, na.rm=TRUE), comp_se=se(complexity), pred=mean(pred_kg, na.rm=TRUE), pred_se=se(pred_kg), pro=mean(po_kg, na.rm=TRUE), pro_se=se(po_kg))%>%left_join(hpqq)
hpqq

write.csv(comp,'final/CSS_means.csv',row.names=FALSE)

summary(aov(po_kg~sp*treat, data=meth))
##predator density differs with plant, treat, interaction (marg)

ggplot(preds, aes(sp, pro, group=treat))+geom_bar(stat='identity',aes(fill=treat),position=position_dodge(.9))+
  geom_errorbar(aes(ymin=pro-pro_se, ymax=pro+pro_se), width=.1, position=position_dodge(.9))

##########
comp
summary(aov(lm(ID~comp,data=comp)))#P=.0435
summary(lm(ID~comp,data=comp))#P=0.40
str(comp)
View(comp)
comp$comp<-ifelse(comp$sp=='SAME', 17.324, comp$comp)
comp$comp<-ifelse(comp$sp=='ARCA', 15.64, comp$comp)

write.csv(comp,'CSS_comp.csv',row.names=FALSE)

ggplot(comp, aes(comp, ID))+
  geom_errorbarh(aes(xmin=comp-comp_se, xmax=comp+comp_se), width=0, alpha=.75, size=.75)+
  geom_errorbar(aes(ymin=ID-ID_var, ymax=ID+ID_var), width=0, alpha=.75, size=.75)+
  geom_smooth(method='lm', se=F, color='darkgrey', lty='solid')+
  geom_hline(aes(yintercept=0), color='grey', lty='dashed')+
  geom_point(size=2.8, stroke=1, fill='black', shape=24)+
  theme_nell()+labs(x='Plant complexity', y='Bird effects\n(LRR control/exclusion)')

##increasing complexity = mre spiders?

cats<-read.csv('cat_region.csv')%>%transform(plant=as.character(plant))%>%
  melt(id.vars=c('sp','plant','location'))%>%
  mutate(value = ifelse(value>=1, 1, 0))%>%
  dcast(sp+plant+variable~location)
View(cats)

prop.test(cats$center, cats$side)

Anova(glm(cbind(center, side, top)~location, data=cats, family='binomial'), type='III')

cat.p<-cats%>%group_by(sp, plant, location)%>%
  summarize(prop = (sum(value)/length(value))/6)

Anova(glm(prop~location, data=cat.p), type='III')

prop.test()
