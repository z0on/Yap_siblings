setwd("~/Dropbox/yap_hetcheck_2020/") 
system("ls *bams*pestPG >ins")
infiles=scan("ins", what="c")

popage=c();pi=c();chr=c()
file=infiles[1]
for(file in infiles){
  tt=read.table(file)
  chr=append(chr,tt[,2])
  tt=tt[,5]/tt[,14] # dividing theta by length = pi
  pi=append(pi,tt)
  popage=append(popage,rep(sub(".bams.+","",file),length(tt)))
}
pop=sub("_.+","",popage)
age=sub(".+_","",popage)

pis=data.frame(cbind(popage,pop,age,chr))
pis$pi=pi  

library(ggplot2)

ggplot(pis,aes(pop,pi,fill=age))+
  geom_violin()+
  theme_bw()

choose="NMP"
gg=ggplot(subset(pis,pop==choose),aes(age,pi))+
  geom_boxplot(outlier.shape=NA)+
  ggtitle(choose)+
  theme_bw()+
  theme(legend.position = "none")+
  ylim(0.0038,0.0058)+
 geom_line(aes(group=chr,color="red"), position = position_dodge(0.1)) +
  geom_point(aes(group=chr),size=2,shape=21, position = position_dodge(0.1))
pdf(paste(choose,"_pi.pdf",sep=""),width=1.7,height=2.5)
plot(gg)
dev.off()


library(lme4) 
library(lmerTest)
ll=lmer(pi~pop+pop:age+(1|chr),pis)
summary(ll)
# Random effects:
#   Groups   Name        Variance  Std.Dev. 
# chr      (Intercept) 1.184e-07 3.441e-04
# Residual             4.346e-09 6.593e-05
# Number of obs: 112, groups:  chr, 14
# 
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)  4.992e-03  9.365e-05  1.384e+01  53.302  < 2e-16 ***
#   popNMP      -2.415e-05  2.492e-05  9.100e+01  -0.969 0.335100    
#   popST       -9.080e-05  2.492e-05  9.100e+01  -3.644 0.000446 ***
#   popWOR      -4.992e-05  2.492e-05  9.100e+01  -2.004 0.048093 *  
#   popGC:ageJ  -5.008e-05  2.492e-05  9.100e+01  -2.010 0.047419 *  
#   popNMP:ageJ -1.734e-04  2.492e-05  9.100e+01  -6.960 5.12e-10 ***
#   popST:ageJ   6.960e-05  2.492e-05  9.100e+01   2.793 0.006363 ** 
#   popWOR:ageJ  4.609e-05  2.492e-05  9.100e+01   1.850 0.067583 .  

library(MCMCglmm)
mm=MCMCglmm(pi~pop+pop:age,random=~chr,data=pis)
HPDinterval(mm$Sol)
# lower         upper
# (Intercept)  4.801758e-03  5.201482e-03
# popNMP      -6.991101e-05  2.257452e-05
# popST       -1.442490e-04 -4.659526e-05
# popWOR      -9.600271e-05 -3.001665e-06
# popGC:ageJ  -9.837307e-05  5.401380e-07
# popNMP:ageJ -2.202942e-04 -1.232545e-04
# popST:ageJ   2.113559e-05  1.209185e-04
# popWOR:ageJ -6.031412e-06  9.127099e-05
# attr(,"Probability")
# [1] 0.95
apply(mm$Sol,2,median)
# (Intercept)        popNMP         popST        popWOR    popGC:ageJ   popNMP:ageJ 
# 5.000101e-03 -2.382581e-05 -9.112022e-05 -4.883477e-05 -5.038738e-05 -1.732624e-04 
# popST:ageJ   popWOR:ageJ 
# 7.037865e-05  4.441474e-05 
         
