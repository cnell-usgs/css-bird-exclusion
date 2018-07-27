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
View(methods)

#arthropod matrix
arths<-read.csv('final/CSS_master_MAT.csv')
arths[is.na(arths)]<-0
View(arths)

##arthropod groupings to determine pred/herbs
orders<-read.csv('final/CSS_feeding.csv')
arth.mat<-arths%>%dplyr::select(ACAR:WASP)##community matrix
sum(rowSums(arth.mat, na.rm=TRUE)) #2789? which plant 1 or 2) are missing?
##2789 arthropods

#################
##compare to 'comp' df - this is the master abundances but has no size data
# want to add in size data for each sample (but some samples missing)
# double check that total abundance = L + S in comp df
# if not, make new column to calcualte difference
# for samples that are not in comp, or named differently - check species, treat, abun to see if misnamed
arth.melt<-arths%>%dplyr::select(Sample:treat)%>%
  melt(id.vars=c('species','Sample','treat'), variable.name='ID', na.rm=TRUE)%>%
  left_join(orders, by='ID')%>%
  dplyr::select(everything(), abun=value)%>%
  transform(ID = as.character(ID), Sample=as.character(Sample))%>%
  mutate(ID_new = ifelse(ID %in% c('ANT','BEE','WASP'), 'HYME', ID))%>%
  group_by(species, Sample, treat, Order = ID_new, feed)%>%
  summarize(abundance = sum(abun, na.rm=TRUE))%>%mutate(orig='arth')%>%
  transform(feed=as.character(feed))%>%ungroup()#fix groupings to match size data
str(arth.melt)
sum(arth.melt$abundance)#2789

##official list of sample and treatments and total abundance
design<-methods%>%select(species, treat, Sample)
str(design)

# this is size data broken down by order - missing some shit
comp<-read.csv('final/CSS_order.csv')%>%
  dplyr::select(species = sp, Sample, ID=Order, value, feed, size)%>%
  transform(ID=as.character(ID))%>%
  mutate(Order = ifelse(ID == 'PSCO', 'PSOC', ID))%>%
  group_by(species, Sample, Order, feed,size)%>%
  summarize(value=sum(value, na.rm=TRUE))
levels(comp$size)<-c('S','L','S')
comp.cast<-comp%>%dcast(species+Sample+Order+feed~size, value.var='value', fun.aggregate=sum)%>%
  mutate(L_S = L+S, orig='comp')%>%transform(feed=as.character(feed), Sample=as.character(Sample))
str(comp.cast)
sum(comp.cast$L_S, na.rm=TRUE)#2091

#left.df<-left_join(arth.melt, comp.cast, by=c('species','Sample','Order'))
#View(left.df)
sum(left.df$L_S, na.rm=TRUE)#1060

full.df<-full_join(arth.melt, comp.cast, by=c('species','Sample','Order','feed'))%>%
  filter(abundance != 0) ##these are only ones that matched to the arth df, missing
View(full.df)

write.csv(sizes, 'final/CSS_arth_sizes.csv', row.names=FALSE)
unique(sizes$Order)
##################
##use sizes to calculate approximate arthopod biomass
sizes<-read.csv('final/CSS_arth_sizes.csv')%>%dplyr::select(-L_S, -orig.x)%>%
  transform(Order=as.character(Order), Sample=as.character(Sample))
arth_b<-read.csv('/Users/colleennell/Dropbox/Projects/CAWR/CAWR_arth/data/cawr_arth_biomass.csv')%>%dplyr::select(Order,a,b,bar)

levels(arth_b$Order)<-c('ACAR','ARAN','ARCH','AUCH','CHIL','COLE','DERM','DIPT','ENTO','HEMI','HETE','HYME','ISOP','LEPI',
                        'MANT','NEUR','OPIL','PHAL','ORTH','PSEU','PSOC','RAPH','SIPH','STER','THYS','THYSA')
arth_b$Order<-as.character(arth_b$Order)
sizes$bm_name<-ifelse(sizes$Order == 'COCC', 'COLE', ifelse(sizes$Order == 'MIRID','HETE', ifelse(sizes$Order == 'ARCH', 'THYS', sizes$Order)))

##calculate biomass - Hodar Equation: W=aBL^b
sized<-sizes%>%left_join(arth_b, by=c('bm_name'='Order'))%>%
  mutate(S_bm = (a*1.5^b)*S, L_bm =(a*4^b)*L, total_bm=S_bm+L_bm)%>%dplyr::select(-bar)
write.csv(sized, 'final/CSS_arthropod_biomass.csv', row.names=FALSE)

##total arthropo biomass on each plant
arth_bm<-sized%>%group_by(species, Sample, treat)%>%
  summarize(sum=sum(total_bm, na.rm=TRUE))

ggplot(arth_bm, aes(interaction(Sample, treat), sum))+geom_text(aes(label=Sample, color=treat))+facet_wrap(~species)

##add to methods df
meth<-methods%>%left_join(arth_bm, by=c('sp'='species','Sample','treat'))%>%select(everything(), arth_biomass=sum)%>%
  mutate(kg_min = plant_kg/vac_min)
write.csv(meth, 'final/CSS_methods_arthropod_biomass_oct27.csv', row.names=FALSE)

##dry biomass vacuumed in 1 min for each plant - arth mg by plant kg
##to compare insect food availability, biomass was calcualted absed ont he collection of insect and plant biomass data 
effort<-meth%>%
  group_by(sp)%>%
  summarize(mean_kg_min = mean(kg_min, na.rm=TRUE), se_kg_min = se(kg_min))#average plant biomass vacummed in min
ggplot(effort, aes(sp, mean_kg_min))+geom_point()+geom_errorbar(aes(ymin=mean_kg_min-se_kg_min, ymax=mean_kg_min+se_kg_min))
#SAAP had most biomass vacuumed per min, ENCA had lowest

##arthropod by sampling effort
dens<-meth%>% #calculate arthopod density
  left_join(effort, by='sp')%>%
  mutate(kg_vac = mean_kg_min*vac_min, ###convert time of vacuuming to plant biomass vacuumed
         arth_kg_vac = abun_comp/kg_vac, #calculate arthropod density,
         bm_kg = arth_biomass/kg_vac) # arth biomass / plant

#arthropod densities
ggplot(dens, aes(interaction(Sample, treat), arth_kg_vac, color=treat))+geom_text(aes(label=Sample))+facet_wrap(~sp, scales='free')

#arthropod biomass
ggplot(dens, aes(interaction(Sample, treat), bm_kg, color=treat))+geom_text(aes(label=Sample))+facet_wrap(~sp, scales='free')

#average density bu species * treatmetn
dens.sum<-dens%>%
  group_by(sp, treat)%>%
  summarize(dens_means = mean(arth_kg_vac, na.rm=TRUE), dens_se = se(arth_kg_vac),
            mg_mean = mean(bm_kg, na.rm=TRUE), mg_se=se(bm_kg))

ggplot(dens.sum, aes(treat, mg_mean))+geom_point()+
  geom_errorbar(aes(ymin=mg_mean-mg_se, ymax=mg_mean+mg_se), width=.2)+
  facet_wrap(~sp)

#SLA - ARCA = 10.5 mm2 mg-1
#leaf size= 51
#########################################
## community composition analyses
###arthropod biomass compostion matrix
library(vegan)
##does bird feeding change the compostion or arthropod communities?

#treat effect on arth biomass?
ord.dens<-sized%>%left_join(meth%>%dplyr::select())


##relative biomass hould use hte plant kg
size.matrix<-sized%>%dcast(species+Sample+treat~Order, value.var='total_bm')
size.matrix[is.na(size.matrix)]<-0
size.mat<-size.matrix%>%select(ACAR,ARAN,AUCH:THYS)
rel.size<-size.mat%>%mutate(total=rowSums(.))%>%mutate_all(funs(./total))%>%select(-total)
size.dist<-vegdist(log(size.mat+1), method='bray')
size.ad<-adonis(size.dist~treat*species, data=size.matrix, permutations=4000)
size.ad
anova(betadisper(size.dist, size.matrix$treat))
#the relative biomass in different arthopod orders differs with treatment, species, and interaction
size.mds<-metaMDS(size.dist, distance='bray',trymax=200)
size.mds
##stress = .18
env.df<-envfit(size.mds, size.matrix%>%dplyr::select(-Sample))
size.df<-data.frame(scores(size.mds), Sample=size.matrix$Sample, sp =size.matrix$species, treat=size.matrix$treat)

#MIRID, HYME, DIPT, THYS, STER, COLE, COCC, AUCH, ARAN
sp.df<-as.data.frame(env.df$factors$centroids)
sp.df$species<-rownames(sp.df)

ord.df<-as.data.frame(env.df$vectors$arrows)
ord.df$order<-rownames(ord.df)
str(ord.df)


ggplot(size.df, aes(NMDS1, NMDS2))+
  theme_nell()+scale_fill_manual(values=c('white','grey'))+
  geom_segment(data=ord.df, aes(x=0, xend=NMDS1, y=0, yend=NMDS2))+
  geom_text(data=sp.df, aes(NMDS1, NMDS2, label=paste(gsub('species','',species)),size=20))#+geom_point(aes(color=sp, fill=treat), shape=21, size=4, stroke=2)

####################################################
##calculate density of herbivores, predators for each plant
feeds<-arth.melt%>%dcast(species+Sample+treat~feed, fun.aggregate=sum)%>%
  mutate(herb_abun = herb+omni, pred_abun = pred+omni)%>%
  left_join(methods)%>%
  mutate(herb_kg = (herb_abun/vac_min)/plant_kg, herb_vol = (herb_abun/vac_min)/plant_vol,
         pred_kg = (pred_abun/vac_min)/plant_kg, pred_vol = (pred_abun/vac_min)/plant_vol,
         arth_kg = (abun_comp/vac_min)/plant_kg, arth_vol = (abun_comp/vac_min)/plant_vol,
         arth_dens = abun_comp/plant_kg, arth_vol_kg = (abun_comp)/plant_vol/plant_kg)

##herbivore densities by plant, treatment
arth.aov<-aov(log(arth_dens+1)~species*treat, data=feeds)
summary(arth.aov)
shapiro.test(resid(arth.aov))

##arthropod densities between treatments
arthdens<-feeds%>%group_by(species, treat)%>%summarize(mean=mean(arth_kg, na.rm=TRUE), se=se(arth_kg), n = length(unique(Sample)))

ggplot(herbdens, aes(x=treat, y= mean))+geom_point()+geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+facet_wrap(~species)

####################################################
##indirect defense by birds
## LRR = log(exclusion/control)
arth.lrr<-arthvolkg%>%melt(id.vas=c('species','treat'))%>%dcast(species~treat+variable)
View(arth.lrr)

library(metafor)
bird.lrr<-escalc('ROM', m1i=C_mean, m2i=T_mean, sd1i=C_sd, sd2i=T_sd, n1i=C_n, n2i=T_n, data=arth.lrr, append=TRUE, var.names=c('LRR','LRR_var'))
bird.ci<-summary(bird.lrr)

##with se
ggplot(bird.ci, aes(x=species, y=LRR))+geom_point()+geom_hline(yintercept=0)+
  geom_errorbar(aes(ymin=LRR-LRR_var, ymax=LRR+LRR_var), width=0)

##with 95% CI
##all but SAME intersect with 0 - i.e. no effect
##but supposed to only use half of these values to calculate LRR?
ggplot(bird.ci, aes(x=species, y=LRR))+geom_point()+geom_hline(yintercept=0)+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub), width=0)

summary(lm(LRR~T_mean, data=bird.lrr))

##should it be a LRR? o rthe total effect of the exclusion treatment - total density consumed
bird.md<-escalc('ROM', m1i=C_mean, m2i=T_mean, sd1i=C_sd, sd2i=T_sd, n1i=C_n, n2i=T_n, data=arth.lrr, append=TRUE)
bird.md.ci<-summary(bird.md)
ggplot(bird.md.ci, aes(x=species, y=yi))+geom_point()+geom_hline(yintercept=0)+
  geom_errorbar(aes(ymin=yi-sei, ymax=yi+sei), width=0)
View(bird.md.ci)

##arthropod density in exclusion
#direct defense
set.seed(56)
#feeds
##randomly pick replicates from exclusiont treatment for DD - 3 or 4 per species
feed.treat<-feeds%>%filter(treat == 'T')%>%
  mutate(use = rep(c('DD','ID'), length.out = length(treat)))
View(feed.treat)
##now calculate density by assigned calculation
feed.treat.m<-feed.treat%>%
  group_by(species, use)%>%
  summarize(vol_dens_mean=mean(arth_vol_kg, na.rm=TRUE), vol_dens_sd=sd(arth_vol_kg, na.rm=TRUE),
            dens_mean=mean(arth_dens, na.rm=TRUE), dens_sd=sd(arth_dens, na.rm=TRUE), n=length(unique(Sample)))%>%
  melt(id.vars=c('species','use'))%>%
  dcast(species~use+variable)%>%left_join(arth.lrr, by='species')
##so instead of using the 'T' category for the bird exlcusion treatment, should use ID/C
lrr.id<-escalc('ROM', m1i=C_mean, m2i=ID_vol_dens_mean, sd1i=C_sd, sd2i=ID_vol_dens_sd, n1i=C_n, n2i=ID_n, data=feed.treat.m, append=TRUE)
lrr.id.ci<-summary(lrr.id)

ggplot(lrr.id.ci, aes(x=log(DD_vol_dens_mean), y=yi))+geom_point()+geom_hline(yintercept=0)+
  geom_errorbar(aes(ymin=yi-sei, ymax=yi+sei), width=0)

lrr.id$mean_diff<-100*(lrr.id$ID_vol_dens_mean-lrr.id$C_mean)/lrr.id$ID_vol_dens_mean
ggplot(lrr.id, aes(log(DD_vol_dens_mean), mean_diff))+geom_point()+
  geom_hline(yintercept=0)

##division seems fine- no differences gerneated from that
summary(aov(arth_vol_kg~use*species, data=feed.treat))

ggplot(feed.treat.m, aes(x=use, y=vol_dens_mean))+
  geom_point()+
  geom_errorbar(aes(ymin=vol_dens_mean-vol_dens_se, ymax=vol_dens_mean+vol_dens_se))+
  facet_wrap(~species)

ggplot(feed.treat.m, aes(x=use, y=dens_mean))+
  geom_point()+
  geom_errorbar(aes(ymin=dens_mean-dens_se, ymax=dens_mean+dens_se))+
  facet_wrap(~species)
