## simulations for Laurell 1994 percentage change and baseline
# x=pre-treatment ppd
# y=post-treatment ppd
library(MASS)
nitn<-100000
nsim<-47# will be 9?
rxy<-0.207 #herb dens in 2 treats = 0.59 =cor(birds$resist, birds$C)
rrx<-0.354 #0.72 = cor(birds$birdfx, birds$resist)
sdx<-2.282 # 0.00967 sd(birds$resist)
sd(birds$birdfx) #0.365

rxy<-0.59
rrx<-0.72

sdx<-sd(birds$resist)
sdy<-sd(birds$C)
meanx<-mean(birds$resist)
meany<-mean(birds$C)

vx<-sdx/meanx
vy<-sdy/meany
vx #0.4436
vy # 0.320
k<-vx/vy
k#1.385
# great id with higher dd vx>vy

# the observed sd of y is 1.406. Under the null hypothesis
# that there is no relation between percentage change and
# baseline, the expected sd of y is 2.282*(3.02/8.45) = 0.816
mua<-c(vx, vy)
sigma<-matrix(c(sdx*sdx, rxy*sdx*sdy, rxy*sdx*sdy, sdy*sdy), nrow=2)
m<-mvrnorm(100000, mua, sigma)
nitn<-10000
nsim<-9
param1<-c(1:nitn)
param2<-c(1:nitn)

for(i in 1:nitn){
  s<-sample(1:nitn,nsim)
  smc<-m[s,]
  param1[i]<-cor(log(smc[,2]/smc[,1]),smc[,2])
  param2[i]<-cor((smc[,1]-smc[,2])/smc[,1],smc[,1])
}
robs<-0.72
p1<-length(param1[param1>=robs])
p1
p1<- p1/nitn
p1 #- likelihood that correlation r is more extreme than observed
## p = 0.0316

p2<-length(param1[param2<=robs])
p2<- 2*p1/nitn
hist(param1,col="grey",main="distribution of correlation
     coefficients",xlab="Correlation coefficients between
     percentage change and baseline",breaks=100)
quantile(-param1, c(.05, .5, .95))
quantile(param2, c(.025, .5, .975))

birds
amod
cmod
emod
raw.cor
sdy<-0.816 #sd(birds$birdfx) #0.365
mua<-c(8.45,3.02)
sigma<-
  matrix(c(sdx*sdx,rxy*sdx*sdy,rxy*sdx*sdy,sdy*sdy),nrow=2)
m<-mvrnorm(100000, mua, sigma)
param1<-c(1:nitn)
param2<-c(1:nitn)
for(i in 1:nitn){
  s<-sample(1:nitn,nsim)
  smc<-m[s,]
  param1[i]<-cor(smc[,2]/smc[,1],smc[,1])
  param2[i]<-cor((smc[,1]-smc[,2])/smc[,1],smc[,1])
}
p1<-length(param1[param1>=-0.354])
p1<- 2*p1/nitn
p2<-length(param1[param2<=0.354])
p2<- 2*p1/nitn
hist(param2,col="grey",main="distribution of correlation
     coefficients",xlab="Correlation coefficients between
     percentage change and baseline",breaks=100)
quantile(param1, c(.025, .5, .975))
quantile(param2, c(.025, .5, .975))
p1
## simulations for CD4 count percentage change and baseline
# x=pre-treatment cd4
# y=post-treatment cd4
library(MASS)
nitn<-100000
nsim<-100
rxy<-0.4835
rrx<-0.7109
sdx<-91.374
# the observed sd of y is 92.44. Under the null hypothesis
# that there is no relation between percentage change and
# baseline, the expected sd of y is
# 91.374*(414.18/258.78) = 146.245
sdy<-146.245
mua<-c(258.78,414.18)
sigma<-
  matrix(c(sdx*sdx,rxy*sdx*sdy,rxy*sdx*sdy,sdy*sdy),nrow=2)
m<-mvrnorm(100000, mua, sigma)
param1<-c(1:nitn)
param2<-c(1:nitn)
for(i in 1:nitn){
  s<-sample(1:nitn,nsim)
  smc<-m[s,]
  param1[i]<-cor(smc[,2]/smc[,1],smc[,1])
  param2[i]<-cor((smc[,1]-smc[,2])/smc[,1],smc[,1])
}
p1<-length(param1[param1>=0.7109])
p1<- 2*p1/nitn
p2<-length(param1[param2<=0.7109])
p2<- 2*p1/nitn
hist(param2,col="grey", xlab="Correlation coefficients between
percentage change and baseline",main="distribution of
correlation coefficients",breaks=100)
quantile(param1, c(.025, .5, .975))
quantile(param2, c(.025, .5, .975)) 