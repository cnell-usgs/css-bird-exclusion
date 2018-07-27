Multivariate- species matrix analyses- PERMANOVA, SIMPER, NMDS

div<-read.csv("mex_all plot.csv")
View(div) #data file with samples as rows and species abundances in columns. can also be applied to plant chemistry analyses and anything that adheres to this format.
library(ggplot2)
library(vegan)

#make data matrix
div.matrix<-vegdist(div[,9:15],method="bray") # can also use method="jaccard" for absense/presence data

#test assumptions- homogeneity of dispersion
div.disp<-betadisper(div.matrix,div$Div)
permutest(div.disp,permutations=999)# output displays mean distance from centroids by groups

#PERMANOVA
div.perm<-adonis(div.matrix~Div,div,permutations=10000,method="bray")
div.perm

#Similarity analysis
div.simper<-simper(div[,9:15],div$Div,permutations=10000,trace=F,parallel=getOption("mc.cores"))
Summary(div.simper) # gives species contributions to overall dissimilarity between groups in rank order
#in your output you can pull 'average' contribution, cumulative sum of dissimilarity 'cusum'

##NMDS plotting
div.nmds<-metaMDS(div.matrix,distance="bray",k=2,trymax=1000)
div.nmds
stressplot(div.nmds) # a good plot has high R2, signifying a good fit between distance matrix and NMDS ordination
#extract MDS scores and put in df
div.plots<-div[,1:2]
div.plots$NMDS1<-div.nmds$points[,1]
div.plots$NMDS2<-div.nmds$points[,2]

#plot in ggplot2
div.fig<-ggplot(div.plots,aes(NMDS1,NMDS2,shape=Div)+geom_point(size=4)+theme_classic()
div.fig
#add 95% CI around treatment groups
div.fig+stat_ellipse(aes(x=NMDS1,y=NMDS2,lty=Div),level=0.95,alpha=0.9,size=1,label=T)
#add species loading vectors
div.sp.vec<-envfit(div.nmds$points,div[,9:15],perm=1000)
div.sp.vec
div.vec<-as.data.frame(div.sp.vec$vectors$arrows*sqrt((div.sp.vec$vectors$r)))
div.vec$species<-rownames(div.sp.vec)
library(grid)
+geom_segment(data=div.vec,aes(x=0,xend=MDS1,y=0,yend=MDS2),arrow=arrow(length=unit(0.5,"cm")),color="grey")+coord_cartesian(xlim=c(-.7,.7),ylim=c(-.6,.6))
div.fig<-ggplot()+geom_point(data=div.plots,aes(NMDS1,NMDS2,shape=Div),size=3)+geom_segment(data=div.vec,aes(x=0,xend=MDS1,y=0,yend=MDS2),arrow=arrow(length=unit(0.5,"cm")),color="grey",inherit_aes=F)+coord_cartesian(xlim=c(-.7,.7),ylim=c(-.6,.6))
div.fig
#add vector labels (there are probably easier and better ways)
div.fig+geom_text(data=div.vec,aes(x=MDS1,y=MDS2,label=c("CA","FL","FR","GR","IN","NE","OM")),size=5,position_dodge=.5)
#to adjust vector lengths (but keep scaled relative to R2 values) to fit in frame
div.fig<-ggplot()+geom_point(data=div.plots,aes(NMDS1,NMDS2,shape=Div),size=3)+geom_segment(data=div.vec,aes(x=0,xend=MDS1*.6,y=0,yend=MDS2*.6),arrow=arrow(length=unit(0.5,"cm")),color="grey")+coord_cartesian(xlim=c(-.7,.7),ylim=c(-.6,.6))+theme_classic()
