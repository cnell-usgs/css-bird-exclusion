###################################  
## Make phylogeny of species

library(dplyr)
library(ggplot2)
setwd("/Users/colleennell/Dropbox/Projects/CSS exclusion/")

## https://github.com/jinyizju/S.PhyloMaker
# seed plant phylogeny 
# Usage
# S.PhyloMaker(splist, tree, nodes, output.splist = T, scenarios = c("S1", "S2", "S3"))
source('https://raw.githubusercontent.com/jinyizju/S.PhyloMaker/master/R_codes%20for%20S.PhyloMaker')
# Arguments
# splist:  data frame which contains a list of seed plant species 
# should include three columns named species, genus, and family
sps<-read.csv('data/2018/CSS_taxa.csv')
View(sps)

# Example
library("phytools") # load the "phytools" package.
example<-read.csv("example.splist.csv",header=T) # read in the example species list.
phylo<-read.tree("PhytoPhylo.tre") # read in the megaphylogeny.
nodes<-read.csv("nodes.csv",header=T) # read in the nodes information of the megaphylogeny.
result<-S.PhyloMaker(spList=example, tree=phylo, nodes=nodes) # run the function S.PhyloMaker.
str(result) # the structure of the ouput of S.PhyloMaker.
par(mfrow=c(1,3),mar=c(0,0,1,0)) # show the phylogenies of the three scenarios.
plot(result$Scenario.1,cex=1.1,main="Scenarion One")
plot(result$Scenario.2,cex=1.1,main="Scenarion Two")
plot(result$Scenario.3,cex=1.1,main="Scenarion Three")

# Qian, H. and Y. Jin. (2016) An updated megaphylogeny of plants, a tool for generating plant phylogenies 
# and an analysis of phylogenetic community structure. Journal of Plant Ecology 9(2): 233–239.

sps<-read.csv('data/2018/CSS_taxa.csv')%>%
  mutate(species = gsub(pattern=' ', replacement = '_', species))

#newick
# from http://phylodiversity.net/phylomatic/
# zanne
css_newick<-'((((((((((((((((((((((((((((((((((((((Artemisia_californica:22.8719,Artemisia_douglasiana:22.8719)artemesia:22.8719,(Isocoma_menziesii:22.8719)isocoma:22.8719,((((((((((((((((((((((((((((((((((((((((((((((((((((Ericameria_palmeri:0.005820):0.000775):0.002797):0.002429):0.008210):0.001060):0.015722):0.004630):0.000887):0.022416):0.051759):0.007346):0.086553):0.028730):0.031924):0.074815):0.050594):0.078045):0.081296):0.156558):0.198377):0.088617):0.150068):0.271884):0.060133):0.482707):0.297078):0.228042):0.244143):1.242861):0.440424):1.219405):0.860231):1.132681):0.254408):0.138646):1.275601):1.030030):0.151178):1.788874):3.584537):4.249170):4.182503):3.988555,(((((((((((((((((((Encelia_californica:1.488245):0.718813):0.133989):0.187249):0.185815):1.491389):2.642231):1.249715):1.243290):1.960068):1.212117):0.941629):0.202329):0.126082):6.001613):1.067463):2.431798):0.509263):0.800646):3.678774):3.968666):2.699318):5.485976):0.436373):1.477749):0.174693):0.483908):1.901311):0.843190)Asteraceae:3.346060):7.104420):16.236000):6.015739):0.292069):4.165166):0.119934):0.640582)Asterales:11.753099):0.900124):0.056880):0.252052):0.788782):3.574670)Campanulidae:1.882210,(((((((((((((((((((((((((((((((((Salvia_apiana:1.585850,((Salvia_mellifera:0.960106):0.003844):0.621900):0.022792):0.226163):0.883191):1.002337):0.057590):0.659959):0.583451):0.904751):1.108534):0.339164):4.014169):4.901322):1.135000):1.380493):3.681908):7.185950):3.211645)Lamiaceae:8.510314):0.621584):0.129732):0.182246):0.432809):0.696615):12.782970):5.153970)Lamiales:11.432543):0.542387):1.629560):0.720309):3.417357):17.813603):4.580345)Lamiidae:1.340875):4.737980)Asteridae:8.908438):0.179371,(((((((Eriogonum_fasciculatum:27.0072)eriogonum:27.0072)Polygonaceae:5.087963):24.304300):2.270850):3.783420)Caryophyllales:25.700064):1.536275)Pentapetalae:1.460940)Superasteridae:0.986915):0.020375,(((((((((((((((((((((((((((((((((((((((((((((Lupinus_albifrons:0.091119):0.012328):0.032723):0.007062):0.035326):0.025867):0.033106):0.003393):0.024307):0.104574):0.206139):0.201899):0.381894):0.309602):0.315696):0.472210):0.213894):0.440529):1.230394):0.338983):0.265617):1.694867):0.929627):2.553105):0.488617):0.911691):7.946996):4.413542):7.091714):4.538196):7.081544):3.656891):7.020039):0.671859):3.836941):7.044528):4.515714)Fabaceae:0.276658):2.378927):0.250678)Fabales:34.946600):6.263220)Fabidae:0.420986):3.855140)Rosidae:0.138911)Superrosidae:1.491856)Gunneridae:17.734960)Eudicotyledoneae:43.727816):7.654960):2.257320):3.495264)Mesangiospermae:23.596800):16.908091):8.728959)Angiospermae:108.965000)Spermatophyta:38.467865):10.085114):4.856365):5.208675):11.765386):17.675263);'
css_tree<-read.newick(text=css_newick)
plot(css_tree)
str(css_tree)

## Phylogenetic independent contrast PIC

# phyloegenetic signal
# 1 = similar 0 = not
# statistical nonindependence among species trait values due to their phylogenetic relatedness
# correct for phylogeny
# PIC is a specialized case of PGLS - phylogenetic generalized least squares






