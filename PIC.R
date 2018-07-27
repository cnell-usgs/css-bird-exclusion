########################
## phylogenetic independent contrasts
########################

# data prep
library(dplyr)
library(reshape2)
library(tidyr)
library(metafor)
# plotting
library(ggplot2)
library(cowplot)
# correlations
library(Hmisc)
library(corrplot)
# PIC
library(caper)
library(phytools)
library(phylosignal)
library(phylobase)

setwd('/Users/colleennell/Dropbox/Projects/CSS exclusion/')
theme_set(theme_minimal())

##############################
## data prep ####
all.plant<-read.csv('data/2018/CSS_plant_data.csv')
hpq<-read.csv('data/2017/CSS_means.csv')%>%dplyr::select(sp, comp, comp_se, hpq, hpq_se)

sp.mean<-all.plant%>%
  group_by(species, treat)%>%
  summarize_at(vars(herb_mg_dens), funs(mean(., na.rm=TRUE), sd(., na.rm=TRUE), se, length))%>%
  melt(id.vars=c('species','treat'))%>%
  dcast(species~treat+variable)%>%
  escalc(measure='ROM', m1i=C_mean, m2i=T_mean, sd1i=C_sd, sd2i=T_sd, n1i=C_length, n2i=T_length, data=.)%>%
  summary()%>%
  left_join(hpq, by=c('species'='sp'))%>%
  mutate(resist=-1*T_mean, birdfx=-1*yi, birdfx_se=-1*sei, resist_se=-1*T_se,)

## prepare dataframe to match with phylogeny
levels(sp.mean$species)<-c("Artemisia_californica", "Artemisia_douglasiana","Encelia_californica","Eriogonum_fasciculatum", "Ericameria_palmeri","Isocoma_menziesii","Lupinus_albifrons","Salvia_apiana","Salvia_mellifera")
sp.lrr<-sp.mean%>%
  transform(sps = as.character(species))%>%
  column_to_rownames(.,var = 'sps')%>%
  separate(col=species, into=c('genus', 'sp_ep'), sep='_', remove=FALSE)

## trait correlations
pairs(sp.lrr%>%dplyr::select('LRR bird effect' = yi, 'Direct defense' = T_mean, C_mean, hpq, comp)%>%transform(hpq=log(1+hpq)))

# pairwise regressions?
hpq_lrr<-summary(aov(lm(yi~log(1+hpq), data=sp.lrr)))
hpq_dd<-summary(aov(lm(T_mean~log(1+hpq), data=sp.lrr)))
dd_lrr<-summary(aov(lm(yi~T_mean, data=sp.lrr)))
comp_lrr<-summary(aov(lm(yi~comp, data=sp.lrr)))

##############################
## PIC ####
## account for phylogenetic relatedness (phylogenetic non-independence) in trait-trait or trait-habitat correlations

## read in tree
# CRAAN -phylogenetic comparative methods -  https://cran.r-project.org/web/views/Phylogenetics.html

## assumptions - Cavalli-Sforza & Edwards 1967, Felsenstein 1973, Cooper et al. 2015
# 1. topology of phylogeny is accurate
# 2. branch lengths are correct
# 3. brownian trait evolution - trait variance increases with time

## test assumptions
# relationships among standardized contrasts and node ehights
# absolute values of standardized contrasts and their sd
# heteroscedasticity in model residuals - Purvis & Harvey 1995

## packages 
library(caper) # caper ref - https://cran.r-project.org/web/packages/caper/vignettes/caper.pdf
library(phytools)
library(ape)

## build phylogeny - zanne et al 2014
# newick
zanne_tree<-'((((((((((((((((((((((((((((((((((((((Artemisia_californica:22.8719,Artemisia_douglasiana:22.8719)artemesia:22.8719,(Isocoma_menziesii:22.8719)isocoma:22.8719,((((((((((((((((((((((((((((((((((((((((((((((((((((Ericameria_palmeri:0.005820):0.000775):0.002797):0.002429):0.008210):0.001060):0.015722):0.004630):0.000887):0.022416):0.051759):0.007346):0.086553):0.028730):0.031924):0.074815):0.050594):0.078045):0.081296):0.156558):0.198377):0.088617):0.150068):0.271884):0.060133):0.482707):0.297078):0.228042):0.244143):1.242861):0.440424):1.219405):0.860231):1.132681):0.254408):0.138646):1.275601):1.030030):0.151178):1.788874):3.584537):4.249170):4.182503):3.988555,(((((((((((((((((((Encelia_californica:1.488245):0.718813):0.133989):0.187249):0.185815):1.491389):2.642231):1.249715):1.243290):1.960068):1.212117):0.941629):0.202329):0.126082):6.001613):1.067463):2.431798):0.509263):0.800646):3.678774):3.968666):2.699318):5.485976):0.436373):1.477749):0.174693):0.483908):1.901311):0.843190)Asteraceae:3.346060):7.104420):16.236000):6.015739):0.292069):4.165166):0.119934):0.640582)Asterales:11.753099):0.900124):0.056880):0.252052):0.788782):3.574670)Campanulidae:1.882210,(((((((((((((((((((((((((((((((((Salvia_apiana:1.585850,((Salvia_mellifera:0.960106):0.003844):0.621900):0.022792):0.226163):0.883191):1.002337):0.057590):0.659959):0.583451):0.904751):1.108534):0.339164):4.014169):4.901322):1.135000):1.380493):3.681908):7.185950):3.211645)Lamiaceae:8.510314):0.621584):0.129732):0.182246):0.432809):0.696615):12.782970):5.153970)Lamiales:11.432543):0.542387):1.629560):0.720309):3.417357):17.813603):4.580345)Lamiidae:1.340875):4.737980)Asteridae:8.908438):0.179371,(((((((Eriogonum_fasciculatum:27.0072)eriogonum:27.0072)Polygonaceae:5.087963):24.304300):2.270850):3.783420)Caryophyllales:25.700064):1.536275)Pentapetalae:1.460940)Superasteridae:0.986915):0.020375,(((((((((((((((((((((((((((((((((((((((((((((Lupinus_albifrons:0.091119):0.012328):0.032723):0.007062):0.035326):0.025867):0.033106):0.003393):0.024307):0.104574):0.206139):0.201899):0.381894):0.309602):0.315696):0.472210):0.213894):0.440529):1.230394):0.338983):0.265617):1.694867):0.929627):2.553105):0.488617):0.911691):7.946996):4.413542):7.091714):4.538196):7.081544):3.656891):7.020039):0.671859):3.836941):7.044528):4.515714)Fabaceae:0.276658):2.378927):0.250678)Fabales:34.946600):6.263220)Fabidae:0.420986):3.855140)Rosidae:0.138911)Superrosidae:1.491856)Gunneridae:17.734960)Eudicotyledoneae:43.727816):7.654960,((((((((((((((((Stipa_pulchra:23.8819)stipa:23.8819)Poaceae:1.789542):4.309137):7.809344):6.449827):13.917930):2.180905):1.372227):4.576783)Poales:33.002000)Commelinidae:11.693860):2.668590):19.828273)Petrosaviidae:0.979827)Nartheciidae:12.924900)Monocotyledoneae:17.016266):2.257320):3.495264)Mesangiospermae:23.596800):16.908091):8.728959)Angiospermae:108.965000)Spermatophyta:38.467865):10.085114):4.856365):5.208675):11.765386):17.675263);'
css_tree<-read.newick(text=zanne_tree)
plot(css_tree)
css_root<-root(css_tree, outgroup='Stipa_pulchra', resolve.root=TRUE)# 0 length edge & too many Nnode (should be length(tip.label-1))

## prepare dataframe to match with phylogeny
levels(sp.lrr$species)<-c("Artemisia_californica", "Artemisia_douglasiana","Encelia_californica","Eriogonum_fasciculatum", "Ericameria_palmeri","Isocoma_menziesii","Lupinus_albifrons","Salvia_apiana","Salvia_mellifera")
rownames(sp.lrr)<-sp.lrr$species
sp.lrr$species<-as.character(sp.lrr$species)
sp.lrr<-sp.lrr%>%separate(col=species, into=c('genus', 'sp_ep'), sep='_', remove=FALSE)


# add dummy row to dataframe for outgroup - Stipa pulchra
data<-sp.lrr
temprow <- matrix(c(rep.int(NA,length(data))),nrow=1,ncol=length(data))
# make it a data.frame and give cols the same names as data
newrow <- data.frame(temprow)
colnames(newrow) <- colnames(data)
# rbind the empty row to data
data <- rbind(data,newrow)
data$species<-ifelse(is.na(data$species),'Stipa_pulchra', paste(data$species))
str(data)

## troubleshooting tree problems
# remove node siwth single descendants
has.singles(css_tree)
tree_collapse<-collapse.singles(css_tree, root.edge=TRUE)
# resolve multifurcations with branchs of 0 legth
css_md<-multi2di(tree_collapse) 
css_md<-makeLabel(css_md)# create labels for empty nodes
## set xero-length branches to 1/100000 total tree length
css_md$edge.length[css_md$edge.length==0]<-max(nodeHeights(css_md))*1e-6

##############################
## correlated trait evolution

# calculate pairwise contrasts at nodes for each trait
# calculate ancestral states
# calculate contrasts among ancestors
# correlate n-1 contrasts for pairs of traits

# make comparative data
data<-data%>%mutate(resist=-1*T_mean, birdfx=-1*yi)
css.cd<-comparative.data(css_md, data, names.col="species", warn.dropped=TRUE, na.omit=FALSE)

#bird effects and direct defense - significant iwth and without PIC
ddid.crunch<-crunch(birdfx~resist, data=css.cd)
summary(ddid.crunch)
caic.table(ddid.crunch)

# host plant quality and bird effect - significant with PIC
hpqyi.crunch<-crunch(birdfx~log(1+hpq), data=css.cd)
summary(hpqyi.crunch) # sig
caic.table(hpqyi.crunch)
summary(lm(birdfx~log(1+hpq), data=data))#not sig

#host platn quality and direct defense - sig with and without PIC
ddhpq.crunch<-crunch(resist~log(1+hpq), data=css.cd)
summary(ddhpq.crunch)
caic.table(ddhpq.crunch)
summary(lm(resist~log(hpq+1), data=data)) # slightly less sig

#host platn quality and comp
ddc.crunch<-crunch(comp~log(1+hpq), data=css.cd)
summary(ddc.crunch)
caic.table(ddc.crunch)
summary(lm(comp~log(hpq+1), data=data)) # slightly less sig

# complexity and indirect defense - becomes insignificant with PIC
compyi.crunch<-crunch(birdfx~comp, data=css.cd)
summary(compyi.crunch) # no sig
caic.table(compyi.crunch)
summary(lm(birdfx~comp, data=data)) #marginal

# complexity and resist - 
compr.crunch<-crunch(resist~comp, data=css.cd)
summary(compr.crunch) # no sig
caic.table(compr.crunch)
summary(lm(resist~comp, data=data)) 

df<-as.data.frame(summary(hpqyi.crunch)$coefficients)%>%mutate(x = 'hpq', y= 'birdfx', mod ='PIC')
df<-as.data.frame(summary(hpqyi.crunch)$coefficients)%>%mutate(x = 'hpq', y= 'birdfx', mod ='PIC')
df<-as.data.frame(summary(hpqyi.crunch)$coefficients)%>%mutate(x = 'hpq', y= 'birdfx', mod ='PIC')
df<-as.data.frame(summary(hpqyi.crunch)$coefficients)%>%mutate(x = 'hpq', y= 'birdfx', mod ='PIC')

str(df)


##############################

## plotting trait correlations
# make figs showing results for lm and corresponding PIC

## plot theme for regressions
css_theme<-list(theme_minimal(),
                theme(axis.line = element_line(color='black'),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      legend.position = 'none', 
                      axis.ticks = element_line(color='black')))

## extract contrast data & test results for figs
# herbivore density and host plant quality
hddy<-ddhpq.crunch$contrast.data$contr$response
hddx<-ddhpq.crunch$contrast.data$contr$explanatory
hdddF<-as.data.frame(cbind(hddx,hddy))
## direct vs indirect defense
ddy<-ddid.crunch$contrast.data$contr$response
ddx<-ddid.crunch$contrast.data$contr$explanatory
dddF<-as.data.frame(cbind(ddx,ddy))
## complexity vs indirect defense
y<-compyi.crunch$contrast.data$contr$response
x<-compyi.crunch$contrast.data$contr$explanatory
dF<-as.data.frame(cbind(x,y))
# complexity and host plant quality
hidy<-hpqyi.crunch$contrast.data$contr$response
hidx<-hpqyi.crunch$contrast.data$contr$explanatory
hiddF<-as.data.frame(cbind(hidx,hidy))


##############################

ddhpq.pic<-ggplot(hdddF, aes(x=hddx,y=hddy))+
  geom_point(size=2.5)+
  labs(y='Herbivore resistance', x='Host Plant Quality')+
  geom_smooth(method='lm', se=FALSE, color='grey40')+
  annotate('text', x=0.006, y=0.0025, label = 'P =  0.001\nR2 = 0.76')+
  scale_y_reverse(limits=c(NA, -0.002))+scale_x_log10()+css_theme

ddhpq.obs<-ggplot(sp.lrr, aes(hpq, T_mean))+
  geom_errorbarh(data=sp.lrr, aes(xmin=hpq-hpq_se, xmax=hpq+hpq_se),color='grey', height=0,size=1)+
  geom_errorbar(aes(ymin=T_mean-T_se, ymax=T_mean+T_se), color='grey',  width=0,size=1)+
  geom_smooth(method='lm', se=F, color='grey40')+
  geom_point()+
  labs(x='Host plant quality', y='Herbivore resistance')+
  
  annotate('text', x=0.14, y=.04, label = 'P = 0.044\nR2 = 0.39')+
  scale_y_reverse()+scale_x_log10()+css_theme


ddid.pic<-ggplot(dddF, aes(x=ddx,y=ddy))+
  geom_point(size=2.5)+
  labs(x='Herbivore resistance', y='LRR Bird effect')+
  geom_smooth(method='lm', se=FALSE, color='grey40')+
  annotate('text', x=0.00275, y=.155, label = 'P = 0.034\nR2 = 0.42')+
  scale_x_reverse()+scale_y_reverse(limits=c(0.18,-.35))+
  css_theme

ddid.obs<-ggplot(sp.lrr, aes(resist, yi))+
  geom_smooth(method='lm', se=F, color='grey40')+
  geom_errorbar(data=sp.lrr, aes(ymin=yi-sei, ymax=yi+sei),color='grey', size=1)+
  geom_errorbarh(aes(xmin=T_mean-T_se, xmax=T_mean+T_se), color='grey', height=0, size=1)+
  geom_point(aes(), size=2.5)+
  labs(x='Herbivore resistance', y='LRR Bird effect')+
  annotate('text', x=0.04, y=.45, label = 'P = 0.028\nR2 = 0.45')+
  scale_x_reverse()+scale_y_reverse()+css_theme
ddid.obs

cid.pic<-ggplot(dF, aes(x=x,y=y))+
  geom_point(size=2.5)+
  labs(x='Structural complexity', y='LRR Bird effect')+
  annotate('text', x=0.45, y=-0.27, label = 'P = 0.18')+scale_y_reverse()+
  css_theme

cid.obs<-ggplot(sp.lrr, aes(comp, yi))+
  geom_smooth(method='lm', se=F, color='grey40', lty='dashed')+
  geom_errorbar(data=sp.lrr, aes(ymin=yi-sei, ymax=yi+sei),color='grey', size=1)+
  geom_errorbarh(aes(xmin=comp-comp_se, xmax=comp+comp_se), color='grey', height=0, size=1)+
  geom_point(aes(), size=2.5)+
  labs(x='Structural complexity', y='LRR Bird effect')+
  annotate('text', x=28, y=-1, label = 'P = 0.089\nR2 = 0.35')+
  scale_y_reverse()+
  css_theme


idhpq.pic<-ggplot(hiddF, aes(x=hidx,y=hidy))+
  geom_point(size=2.5)+
  labs(y='LRR Bird effect', x='Host Plant Quality')+
  geom_smooth(method='lm', se=FALSE, color='grey40')+
  annotate('text', x=0.008, y=-0.25, label = 'P =  0.048\nR2 = 0.37')+
  scale_y_reverse()+scale_x_log10()+css_theme

idhpq.obs<-ggplot(sp.lrr, aes(hpq, yi))+
  geom_errorbarh(data=sp.lrr, aes(xmin=hpq-hpq_se, xmax=hpq+hpq_se),color='grey', height=0,size=1)+
  geom_errorbar(data=sp.lrr, aes(ymin=yi-sei, ymax=yi+sei),color='grey', size=1, width=0)+
  geom_point()+
  labs(x='Host plant quality', y='LRR Bird effect')+
  annotate('text', x=0.16, y=.45, label = 'P = 0.73')+
  scale_y_reverse()+scale_x_log10()+css_theme


ddhpq.dub<-plot_grid(ddhpq.obs, ddhpq.pic, labels=c('observed', 'PIC'), label_size=12, label_x=c(.08,.19))
ddid.dub<-plot_grid(ddid.obs, ddid.pic, labels=c('observed', 'PIC'), label_size=12, label_x=c(.08,.16))
cid.dub<-plot_grid(cid.obs, cid.pic, labels=c('observed', 'PIC'), label_size=12, label_x=c(.08,.16))
idhpq.dub<-plot_grid(idhpq.obs, idhpq.pic, labels=c('observed', 'PIC'), label_size=12, label_x=c(.08,.19))
alldubs<-plot_grid(ddhpq.dub,idhpq.dub,ddid.dub, cid.dub, nrow=4)
alldubs


## this is written such that....
# high herbivore resistance is a low density in exclusion
# negative bird effect /high is stornger effect of birds i.e. greater removal effect
##############################
## pplot contrasts at nodes of phylogeny
plot(css_md)
nodelabels(round(hpqyi.crunch$contrast.data$contr$explanatory [,1], 3), adj = c(0, -0.5), frame="n") ##HPQ


# sp.lrr on species traits can be inclede in phylo regressiont o estimate the relative contribution of traits vs phylogeny ni explaining interactions
directd<-sp.lrr$T_mean
indirectd<-sp.lrr$yi
hostpq<-log(1+sp.lrr$hpq)
complex<-sp.lrr$comp

names(directd)<-row.names(sp.lrr)
names(indirectd)<-row.names(sp.lrr)
names(hostpq)<-row.names(sp.lrr)
names(complex)<-row.names(sp.lrr)

# remove tree tips for missing sp.lrr
sp.out<-css_md$tip.label[!(css_md$tip.label %in% sp.lrr$species)]
css_drop<-drop.tip(css_md, sp.out)
rownames(sp.lrr)<-sp.lrr$species
head(sp.lrr)

# calculate contrasts for each trait - scaled using expected variances
dd.pic<-pic(directd, css_drop, scaled=TRUE, var.contrasts=TRUE)
id.pic<-pic(indirectd, css_drop, scaled=TRUE, var.contrasts=TRUE)
hpq.pic<-pic(hostpq, css_drop, scaled=TRUE, var.contrasts=TRUE)
comp.pic<-pic(complex, css_drop, scaled=TRUE, var.contrasts=TRUE)

# multiple traits at once
View(sp.lrr)

trait.df<-as.data.frame(sp.lrr)%>%dplyr::select(birdfx=yi*-1, birdfx_se=-1*sei, resist = -1*T_mean,resist_se=-1*T_se,  hpq,hpq_se, comp, comp_se, species)
str(as.data.frame(sp.lrr))
sp.lrr$yi
str(trait.df)
cont.css <- as.data.frame(apply(trait.df, 2, pic, css_drop))
View(cont.css)

pairs(cont.css)

plot_corr<-function(xvar, yvar, )

## plot contrasts on phylogeny
plot(css_drop)
nodelabels(round(id.pic, 2))

trait_trait<-function(xvar, yvar, data){
  cor.out<-cor.test(x=xvar, y=yvar, method='pearson')
  
  xlab<-ifelse(xvar == hpq, 'Host plant quality', ifelse(xvar == resist, 'Herbivore resistance', ifelse(xvar == comp, 'Structural complexity', ifelse(yvar == birdfx, 'LRR bird effect', 'error'))))
  ylab<-ifelse(yvar == hpq, 'Host plant quality', ifelse(yvar == resist, 'Herbivore resistance', ifelse(yvar == comp, 'Structural complexity', ifelse(yvar == birdfx, 'LRR bird effect', 'error'))))
  
  cor.fig<-ggplot(data, aes(xvar, yvar))+
    geom_point()+
    labs(x=xlab, y=ylab)
  if(cor.out$p.value <= 0.05){
    linet<-'solid'
    cor.fig<-cor.fig+geom_smooth(method='lm', se=FALSE, color='grey40'. lty=linet)
  }else if(cor.out$p.value <= 0.10){
    linet<-'dashed'
    cor.fig<-cor.fig+geom_smooth(method='lm', se=FALSE, color='grey40'. lty=linet)
  }else{
    cor.fig<-cor.fig
  }
  return(list(corr=cor.out, fig = cor.fig))
}


dd_id_pic<-lm(birdfx~resist-1, data=cont.css)
summary.lm(dd_id_pic)

# plotting data

trait.df<-as.data.frame(sp.lrr)%>%dplyr::select(birdfx=yi*-1, birdfx_se=-1*sei, resist = -1*T_mean,resist_se=-1*T_se,  hpq,hpq_se, comp, comp_se, species)
str(as.data.frame(sp.lrr))
sp.lrr$yi
str(trait.df)
cont.css <- as.data.frame(apply(trait.df, 2, pic, css_drop))

# PIC
ggplot(cont.css, aes(resist,birdfx))+
  geom_point()+
  geom_smooth(method='lm', se=FALSE, color='grey40')+
  labs(x='Herbivore resistance', y='LRR bird effect')+
  css_theme

# observed

resist_bird<-ggplot(data, aes(resist,birdfx))+
  geom_hline(yintercept=0, lty='dashed')+
  geom_smooth(method='lm', se=FALSE, color='grey40')+
  labs(x='Herbivore resistance', y='LRR bird effect')+
  geom_errorbar(aes(ymin=birdfx-birdfx_se, ymax=birdfx+birdfx_se), width=0, color='grey')+
  geom_errorbarh(aes(xmin=resist-resist_se, xmax=resist+resist_se), height=0, color='grey')+
  geom_point()+
  css_theme
resist_bird
hpq_bird<-ggplot(data, aes(hpq,birdfx))+
  geom_hline(yintercept=0, lty='dashed')+
  labs(x='Host plant quality', y='LRR bird effect')+
  geom_errorbar(aes(ymin=birdfx-birdfx_se, ymax=birdfx+birdfx_se), width=0, color='grey')+
  geom_errorbarh(aes(xmin=hpq-hpq_se, xmax=hpq+hpq_se), height=0, color='grey')+
  geom_point()+
  scale_x_log10()+
  css_theme
hpq_bird
comp_bird<-ggplot(data, aes(comp,birdfx))+
  geom_hline(yintercept=0, lty='dashed')+
  geom_smooth(method='lm', se=FALSE, color='grey40', lty='dashed')+
  geom_errorbarh(aes(xmin=comp-comp_se, xmax=comp+comp_se), height=0, color='grey')+
  geom_errorbar(aes(ymin=birdfx-birdfx_se, ymax=birdfx+birdfx_se), width=0, color='grey')+
  labs(x='Structural complexity', y='LRR bird effect')+
  geom_point()+
  css_theme
comp_bird
comp_resist<-ggplot(data, aes(comp,resist))+
  geom_smooth(method='lm', se=FALSE, color='grey40')+
  labs(x='Structural complexity', y='Herbiovre resistance')+
  geom_errorbarh(aes(xmin=comp-comp_se, xmax=comp+comp_se), height=0, color='grey')+
  geom_errorbar(aes(ymin=resist-resist_se, ymax=resist+resist_se), width=0, color='grey')+
  geom_point()+
  css_theme
comp_resist
hpq_resist<-ggplot(data, aes(hpq,resist))+
  geom_smooth(method='lm', se=FALSE, color='grey40')+
  labs(x='Host plant quality', y='Herbivore resistance')+
  geom_errorbar(aes(ymin=resist-resist_se, ymax=resist+resist_se), width=0, color='grey')+
  geom_errorbarh(aes(xmin=hpq-hpq_se, xmax=hpq+hpq_se), wheight=0, color='grey')+
  geom_point()+
  scale_x_log10()+
  css_theme
hpq_resist
hpq_comp<-ggplot(data, aes(hpq,comp))+
  labs(x='Host plant quality', y='Structural complexity')+
  geom_errorbar(aes(ymin=comp-comp_se, ymax=comp+comp_se), width=0, color='grey')+
  geom_errorbarh(aes(xmin=hpq-hpq_se, xmax=hpq+hpq_se), height=0, color='grey')+
  geom_point()+
  scale_x_log10()+
  css_theme
hpq_comp

plot_grid(hpq_resist, comp_resist, hpq_bird, comp_bird, nrow=2, ncol=2, align='vh')

tb<-caic.table(ddhpq.crunch)
plot(resist~`log(1 + hpq)`, data=tb)

plot_grid(hpq_resist, comp_resist, nrow=1)


corr<-rcorr(as.matrix(cont.css%>%dplyr::select(birdfx, resist, comp, hpq)))
corrplot(corr$r, method='color', type='upper', p.mat=corr$P, sig.level =0.05, insig='blank')
corr$P

str(cont.css)
hpq_dd_pic<-lm(dd.pic~hpq.pic-1)
summary.lm(hpq_dd_pic)

dd_id_pic<<-lm(yi~T_mean-1, data=cont.css)
summary(dd_id_pic<)

hpq_dd_pic<-lm(T_mean~hpq-1, data=cont.css)
summary(hpq_dd_pic)
# the -1 specifies that the regressionis through the origin (intercept =0) - Garland et al 1992
# non significant with PIC - observed cofrelation was artefact of evoltionary histories
# 
?barplot.phylo4d

tree4d<-as(css_drop, 'phylo4d')
class(tree4d)
barplot(tree4d)

xx<-contMap(css.cd$phy,css.cd$data$birdfx,plot=FALSE )


plotTree(css_drop)
plotTree.barplot(css_drop, css.cd$data$resist)
??trait.plot



