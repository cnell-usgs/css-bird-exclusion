## css final for manuscript
library(tidyverse)
library(ape)
library(phytools)
library(caper)
library(reshape2)
library(cowplot)
library(broom)

setwd('/Users/colleennell/Dropbox/Projects/CSS exclusion/')
source('css_funs.R')
########################
## read in data ####
#plants<-read.csv('data/2018/CSS_plants.csv')%>%group_by(species)%>%dplyr::summarize(mean=mean(complexity/diameter))
plants<-read.csv('data/2018/CSS_comps.csv')%>%group_by(species)%>%dplyr::summarize(complex=mean(complexity/(3*diameter)), complex_se=se(complexity/(3*diameter)))
plants

## S.PhyloMaker function to generate phylogeny for seed plants
source("https://raw.githubusercontent.com/jinyizju/S.PhyloMaker/master/R_codes%20for%20S.PhyloMaker")
# Citation: Qian, H. and Y. Jin. (2016) An updated megaphylogeny of plants, a tool for generating plant phylogenies and an analysis of phylogenetic community structure. Journal of Plant Ecology 9(2): 233–239.
# uses PhytoPhylo species-level megaphylogeny as a backbone (Zanne et al 2014)
phylo<-read.tree("data/trees/QianJin_2016.txt") # megaphylogeny from Qian & Jim 2016 "PhytoPhylo"
nodes<-read.table('https://raw.githubusercontent.com/jinyizju/S.PhyloMaker/master/nodes', fill=TRUE, header=TRUE)# nodes for phylogeny

sp.list<-read.csv('data/2018/css_taxa.csv')%>%
  dplyr::select(species, genus, family)

# artemisia species are not in megaphy, use representative taxa from subgenera for divergence
new.sp.list<-sp.list%>%
  mutate(species=ifelse(species == 'Artemisia californica', 'Artemisia tridentata', 
                        ifelse(species == 'Artemisia douglasiana', 'Artemisia ludoviciana', paste(species))))
result<-S.PhyloMaker(spList=new.sp.list, tree=phylo, nodes=nodes)# prune megaphy to species list
phy<-result$Scenario.1%>%makeLabel()
#######################
# ratio of predators to herbivores by plant species
oplant<-read.csv('data/2018/CSS_plant_data.csv')
ph.df<-oplant%>%group_by(species, treat)%>%summarize(mean=mean(pred_herb, na.rm=TRUE), se=se(pred_herb))
ggplot(ph.df, aes(treat, mean))+geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+geom_point()+facet_wrap(~species)
mod<-aov(pred_herb~species*treat, data=oplant)

Anova(mod, type='III')
TukeyHSD(mod, 'treat')
#######################
all.plant<-read.csv('data/2018/CSS_plant_data.csv')%>%filter(treat=='T')%>%
  group_by(species)%>%summarize(pred=mean(pred_mg_dens, na.rm=TRUE), se(pred_mg_dens), predherb=mean(pred_herb, na.rm=TRUE), predherb_se=se(pred_herb))
all.plant## correct

## species means for all traits
sp.lrr<-read.csv('data/2018/css_trait_means.csv')%>%
  mutate(species=ifelse(X == 'Artemisia_californica', 'Artemisia_tridentata',
                        ifelse(X == 'Artemisia_douglasiana', 'Artemisia_ludoviciana', paste(X))))%>%
  dplyr::select(T_mean,species,genus, sp_ep, contains('birdfx'), contains('resist'), contains('hpq'), contains('comp'), contains('ph'), lrr_pred, lrr_pred_se)%>%
  mutate(herbs=log(1+T_mean), hpq_log=log(hpq+1), hpq_log_se=log(hpq_se+1), resist_log=log(resist+1), resist=resist*-1)
rownames(sp.lrr)<-sp.lrr$species
#write.csv(sp.lrr, 'CSS_sp_lrr_all.csv')
#write.csv(sp.lrr, 'R/css-bird-exclusion/data/ms/css_trait_means_figs_check.csv')



sp.lrr
exp(mean(sp.lrr$birdfx))
mod<-lm(ph~complex, data=sp.lrr)
summary(mod)
summary(aov(lm(lrr_pred~hpq_log, data=sp.lrr)))
shapiro.test(resid(mod))
############
##plant-level
planty<-read.csv('data/2018/CSS_plant_data_var.csv')
planty
all.plant
all.plant<-read.csv('data/2018/CSS_plant_data.csv')%>%filter(treat=='C')%>%
  group_by(species,treat)%>%summarize(pred=mean(pred_mg_dens, na.rm=TRUE), se(pred_mg_dens), predherb=mean(pred_herb, na.rm=TRUE), predherb_se=se(pred_herb))
all.plant## correct

## mean by treatment -overall
all.plant%>%group_by(treat)%>%dplyr::summarize(mean=mean(herb_mg_dens, na.rm=TRUE), se=se(herb_mg_dens))
## percent difference

## mean density for each species in control and exclusion
sp.mean<-all.plant%>%group_by(species, treat)%>%dplyr::summarize(mean=mean(herb_mg_dens, na.rm=TRUE), se=se(herb_mg_dens))
sp.mean

# forget treatment
not<-all.plant%>%group_by(species)%>%dplyr::summarize(mean=mean(herb_mg_dens, na.rm=TRUE), se=se(herb_mg_dens))
not
max(not$mean)
min(not$mean)

sp.mean%>%dcast(species~treat, value.var='mean')%>%mutate(times=`T`/C, diff = `T`-C, redux = diff/C, LRR = log(`T`/C))
# redux is the difference across all plants
## use plantmeans and take grand mean - variabl
sp.mean%>%dcast(species~treat, value.var='mean')%>%mutate(times=`T`/C, diff = `T`-C, redux = diff/C)%>%
  summarize_if(is.numeric, funs(mean, se))
log(.483)

phy
sp.lrr
#########################
# characterizing indirect defense
css.phy<-comparative.data(phy, sp.lrr, names.col="species", warn.dropped=TRUE, na.omit=TRUE, vcv=TRUE, vcv.dim=3)
phy.05<-compute.brlen(css.phy$phy, power=0.5) # branch length adjustment for PIC
css.cd<-comparative.data(phy.05, sp.lrr, names.col="species", warn.dropped=TRUE, na.omit=TRUE, vcv=TRUE, vcv.dim=3)

amod<-crunch(birdfx~resist, css.cd,names='species', stand.contr=TRUE, equal.branch.length=FALSE)#
caic.diagnostics(amod, plot=TRUE, plot.signif=TRUE) # OK

bmod<-crunch(birdfx~hpq_log, css.cd, stand.contr=TRUE, equal.branch.length=FALSE)#
caic.diagnostics(bmod, plot=TRUE, plot.signif=TRUE) # AGE

cmod<-crunch(resist~hpq_log, css.cd, stand.contr=TRUE, equal.branch.length=FALSE)#
caic.diagnostics(cmod, plot=TRUE, plot.signif=TRUE) # AGE

dmod<-crunch(birdfx~complex, css.cd, stand.contr=TRUE, equal.branch.length=FALSE)#
caic.diagnostics(dmod, plot=TRUE, plot.signif=TRUE) # OK

emod<-crunch(resist~complex, css.cd, stand.contr=TRUE, equal.branch.length=FALSE)#
caic.diagnostics(emod, plot=TRUE, plot.signif=TRUE) # OK

fmod<-crunch(hpq~complex, css.cd, stand.contr=TRUE, equal.branch.length=FALSE)#
caic.diagnostics(fmod, plot=TRUE, plot.signif=TRUE) # OK

gmod<-crunch(ph~complex, css.cd, stand.contr=TRUE, equal.branch.length=FALSE)#
gmod
caic.diagnostics(gmod, plot=TRUE, plot.signif=TRUE) # OK

sp.lrr
raw.cor<-rbind(tidy(cor.test(sp.lrr$birdfx, sp.lrr$resist))%>%mutate(yvar = 'birdfx', xvar = 'resist'),
               tidy(cor.test(sp.lrr$birdfx, sp.lrr$hpq_log))%>%mutate(yvar = 'birdfx', xvar = 'hpq_log'),
               tidy(cor.test(sp.lrr$resist, sp.lrr$hpq_log))%>%mutate(yvar = 'resist', xvar = 'hpq_log'),
               tidy(cor.test(sp.lrr$birdfx, sp.lrr$complex))%>%mutate(yvar = 'birdfx', xvar = 'comp'),
               tidy(cor.test(sp.lrr$resist, sp.lrr$complex))%>%mutate(yvar = 'resist', xvar = 'comp'),
               tidy(cor.test(sp.lrr$hpq_log, sp.lrr$complex))%>%mutate(yvar = 'hpq_log', xvar = 'comp'),
               tidy(cor.test(sp.lrr$ph, sp.lrr$hpq_log))%>%mutate(yvar = 'predherb', xvar = 'hpq_log'),
               tidy(cor.test(sp.lrr$ph, sp.lrr$complex))%>%mutate(yvar = 'predherb', xvar = 'comp'))%>%
  mutate(model='raw', p.one.tail = p.value/2)
raw.cor

atab<-caic.table(amod)
btab<-caic.table(bmod)
ctab<-caic.table(cmod)
dtab<-caic.table(dmod)
etab<-caic.table(emod)
ftab<-caic.table(fmod)
gtab<-caic.table(gmod)

## using grafen branch length transformation, rho=.5 (sqrt)
pic.cor<-rbind(tidy(summary(amod))%>%mutate(yvar = 'birdfx'),
               tidy(summary(bmod))%>%mutate(yvar = 'bird_fx'),
               tidy(summary(cmod))%>%mutate(yvar = 'resist'),
               tidy(summary(dmod))%>%mutate(yvar = 'birdfx'),
               tidy(summary(emod))%>%mutate(yvar = 'resist'),
               tidy(summary(fmod))%>%mutate(yvar = 'hpq_log'),
               tidy(summary(gmod))%>%mutate(yvar = 'ph'))%>%
  mutate(model='pic', p.one.tail = p.value/2)
pic.cor

#######################
birds<-css.cd$data%>%mutate(species=paste0(genus,'\n', sp_ep))

## birdfx
birds$exp_bfx<-exp(birds$birdfx)
birds$C<-birds$resist/birds$exp_bfx ## find C from resist and exp(LRR)
birds$diff<-birds$resist-birds$C # difference in density between treatments
birds$perc_diff<-birds$diff/birds$resist # percent difference
birds

# plotting theme
css_theme<-list(theme_minimal(),
                theme(axis.line = element_line(color='black'),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      legend.position = 'none', 
                      axis.ticks = element_line(color='black'),
                      axis.title = element_text(size=12),
                      axis.text = element_text(size=12)))
# birdfx by shrub
sps.id<-ggplot(birds, aes(reorder(species, birdfx), birdfx))+
  geom_errorbar(aes(ymin=birdfx-birdfx_se, ymax=birdfx+birdfx_se), width=0, color='darkgrey', size=1.2)+
  geom_point(shape=25, size=2.5, fill='white', stroke=1)+
  geom_hline(yintercept=0, lty='dotted')+
  css_theme+
  theme(axis.text.x=element_text(angle=90, size=14, hjust=.95, vjust=.5, face='italic'))+
  labs(x='', y='Indirect defense from birds')
sps.id

#resistnace by shrub
sps.dd<-ggplot(birds, aes(reorder(species, resist), resist))+
  geom_errorbar(aes(ymin=resist-resist_se, ymax=resist+resist_se), width=.15)+
  geom_point()+
  css_theme+scale_y_reverse()+
  theme(axis.text.x=element_text(angle=90, hjust=.95, vjust=.5, face='italic'))+
  labs(x='', y='Herbivore resistance')

#complexity by sps
sps.comp<-ggplot(birds, aes(reorder(species, complex), complex))+
  geom_errorbar(aes(ymin=complex-complex_se, ymax=complex+complex_se), width=.15)+
  geom_point()+
  css_theme+
  theme(axis.text.x=element_text(angle=90, hjust=.95, vjust=.5, face='italic'))+
  labs(x='Shrub species', y='Structural complexity')



plot_grid(sps.id, sps.dd, sps.comp, ncol=1, nrow=3, labels=c('a','b','c'), label_size=10)
birds


##########################
## mean herbivore density in control and exclusiont treatments
# exclusion
mean(birds$T_mean) # this is the same as resistance
# control
mean(birds$herbs)
## take mean of all plants - not species

## indirect defense reduced herbivore densities by
birds
se(birds$birdfx)

sp.mean%>%dcast(species~treat, value.var='mean')%>%mutate(times=`T`/C, diff = `T`-C, redux = diff/C, LRR = log(`T`/C))

# redux is the difference across all plants
## use plantmeans and take grand mean - variabl
sp.mean%>%dcast(species~treat, value.var='mean')%>%mutate(times=`T`/C, diff = `T`-C, redux = diff/C)%>%
  summarize_if(is.numeric, funs(mean, se))
log(.483)

#########################
## figures - raw correlations


sps.id<-ggplot(birds, aes(reorder(species, birdfx), birdfx))+
  geom_errorbar(aes(ymin=birdfx-birdfx_se, ymax=birdfx+birdfx_se), width=0, color='darkgrey', size=1.1)+
  geom_point(shape=25, size=2.5, fill='white', stroke=1)+
  geom_hline(yintercept=0, lty='dotted')+
  css_theme+
  theme(axis.text.x=element_text(angle=90, size=10, hjust=.95, vjust=.5, face='italic'))+
  labs(x='', y='Indirect defense from birds')
sps.id


# Figure 1a - direct vs indirect defense
fig1a<-ggplot(css.cd$data, aes(resist,birdfx))+
  geom_hline(yintercept=0, lty='dotted')+
  geom_smooth(method='lm', se=FALSE, color='grey40', lty='solid')+
  labs(x=' ', y=' ')+
  geom_errorbar(aes(ymin=birdfx-birdfx_se, ymax=birdfx+birdfx_se), width=0, color='darkgrey', size=1, alpha=.7)+
  geom_errorbarh(aes(xmin=resist-resist_se, xmax=resist+resist_se), height=0, color='darkgrey', size=1, alpha=.7)+
  geom_point(size=2.5, shape=21, fill='white', stroke=1)+scale_x_reverse()+
  css_theme+xlim(0.048,0.005)+theme(axis.text.y=element_blank())
plot_grid(sps.id, fig1a, align='h', rel_widths=c(.8,1))

# fig 1b = complexity & resistance
fig1b<-ggplot(css.cd$data, aes(complex,resist))+
  labs(x='Structural complexity', y='Herbivore resistance')+
  geom_smooth(method='lm', se=FALSE, color='grey40')+
  geom_errorbar(aes(ymin=resist-resist_se, ymax=resist+resist_se), width=0, color='grey', size=1, alpha=.7)+
  geom_errorbarh(aes(xmin=complex-complex_se, xmax=complex+complex_se), height=0, color='grey', size=1, alpha=.7)+
  geom_point(size=2, shape=21, fill='white', stroke=1)+scale_y_reverse()+
  css_theme
fig1b

# fig 1c - complexlexity x ID
fig1c<-ggplot(css.cd$data, aes(complex,birdfx))+
  labs(x='Structural complexity', y='Indirect defense from birds')+
  geom_smooth(method='lm', se=FALSE, color='grey40', lty='solid')+
  geom_errorbar(aes(ymin=birdfx-birdfx_se, ymax=birdfx+birdfx_se), width=0, color='grey', size=1, alpha=.7)+
  geom_errorbarh(aes(xmin=complex-complex_se, xmax=complex+complex_se), height=0, color='grey', size=1, alpha=.7)+
  geom_point(size=2.5, shape=21, fill='white', stroke=1)+geom_hline(yintercept=0, lty='dotted')+
  css_theme
fig1c

plot

#complexity and pred:herb
fig1d<-ggplot(sp.lrr, aes(complex, ph))+
  labs(x='Structural complexity', y='Predator:herbivore')+
  geom_smooth(method='lm', se=FALSE, color='grey40', lty='solid')+
  geom_errorbar(aes(ymin=ph-ph_se, ymax=ph+ph_se), width=0, color='grey', size=1, alpha=.7)+
  geom_errorbarh(aes(xmin=complex-complex_se, xmax=complex+complex_se), height=0, color='grey', size=1, alpha=.7)+
  geom_point(size=2, shape=21, fill='white', stroke=1)+
  css_theme
fig1d
gmod
summary(lm(css.cd$data$ph~css.cd$data$resist))
pic.ph<-ggplot(gtab, aes(complex, ph))+
  geom_smooth(method='lm', se=FALSE, color='grey40', formula=y~x-1)+
  geom_point(size=2)+
  labs(x='Structural complexity',y='Predator:herbivore')+
  css_theme
pic.ph
gtab
plot_grid(fig1d, pic.ph, ncol=2, nrow=1, labels=c('a','b'), label_size=10)
pic.cor
raw.cor
#plot trio stacked
## what about statistical significance?
plot_grid(fig1a, fig1b, fig1c, ncol=1, nrow=3, labels=c('a','b','c'), label_size=10)
plot_grid(fig1b, fig1c,fig1d, ncol=1, nrow=3, labels=c('a','b','c'), label_size=10)
fig1a

## Supplement - PIC correlations
pic.birdr<-ggplot(atab, aes(resist, birdfx))+
  geom_smooth(method='lm', se=FALSE, color='grey40', formula=y~x-1)+
  geom_point(size=2)+
  labs(x='Herbivore resistance',y='Indirect defense from birds')+
  css_theme+scale_x_reverse()

pic.hpqr<-ggplot(ctab, aes(hpq_log, resist))+
  geom_smooth(method='lm', se=FALSE, color='grey40', formula=y~x-1)+
  geom_point(size=2)+scale_y_reverse()+
  labs(x='Host plant quality',y='Herbivore resistance')+
  css_theme

pic.complexr<-ggplot(etab, aes(complex, resist))+
  geom_smooth(method='lm', se=FALSE, color='grey40', formula=y~x-1)+
  geom_point(size=2)+
  labs(x='Structural complexity',y='Herbivore resistance')+
  css_theme+scale_y_reverse()

pic.birdc<-ggplot(dtab, aes(complex, birdfx))+
  geom_smooth(method='lm', se=FALSE, color='grey40', formula=y~x-1)+
  geom_point(size=2)+
  labs(x='Structural complexity',y='Indirect defense from birds')+
  css_theme

pic.birdr
pic.complexr
pic.hpqr
plot_grid(pic.birdr, pic.complexr, pic.birdc, ncol=1, nrow=3, labels=c('a','b','c'), label_size=10)
pic.cor
#########
## correlation coefficients

# direct and indirect
amod #pull rsq for correlation coefficient
## PIC r =
sqrt(.3453)
# raw
summary(lm(birdfx~resist,css.cd$data))


############
# fig 2 - ab = hpq vs resist
fig3.raw<-ggplot(css.cd$data, aes(hpq,resist))+
  labs(x='Host plant quality (mg)', y='Herbivore resistance')+
  geom_smooth(method='lm', se=FALSE, color='grey40')+
  geom_errorbar(aes(ymin=resist-resist_se, ymax=resist+resist_se), width=0, color='grey', size=1, alpha=.7)+
  geom_errorbarh(aes(xmin=hpq-hpq_se, xmax=hpq+hpq_se), height=0, color='grey', size=1, alpha=.7)+
  geom_point(size=2)+scale_y_reverse()+
  scale_x_log10()+
  css_theme

plot_grid(fig3.raw, pic.hpqr, nrow=1, ncol=2, labels=c('a','b'))

ctab
summary(cmod)
sqrt(.6655)
summary(lm(birdfx~resist,css.cd$data))
cor.test(css.cd$data$birdfx, css.cd$data$complex)


emod
sqrt(.5431)


