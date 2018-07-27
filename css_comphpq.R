###quick look at plant complexity and hpq from CSS experiment
#complexity measures amy be inconsistent between years?
###species that were rpeviously measureed were remeasured (mostly for species that could be relocated)
###so this is real dirty

setwd("/Users/colleennell/Documents/R/CSS/data")
css<-read.csv("CSS_complexity.csv")
View(css)

css$hpq.log<-log(css$hpq)
css$hpq.log
css[css=="-Inf"]<-0
lm1<-lm(hpq~comp2,data=css)
summary(lm1)
library(ggplot2)
comp.plot<-ggplot(css,aes(x=comp2,y=hpq))+geom_point(size=4)+geom_smooth(method="lm",se=F)+theme_minimal()+
  geom_errorbar(aes(ymin=hpq-hpq_se,ymax=hpq+hpq_se))+geom_errorbarh(aes(xmin=comp2-comp_se2,xmax=comp2+comp_se2))+
  labs(x="Complexity (Intersections/length)",y="Host Plant Quality (caterpillar wt gain)")+
  geom_text(aes(label=Plant.ID),hjust=-.30,vjust=0.3)
comp.plot

y<-x[!is.na(x)]