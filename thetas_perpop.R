setwd("~/Dropbox/yap_hetcheck_2020/Yap_siblings") 
system("ls *popbams.pestPG >ins")
infiles=scan("ins", what="c")

popage=c();pi=c();chr=c()
file=infiles[6]
for(file in infiles){
  tt=read.table(file)
  chr=append(chr,tt[,2])
  tt=tt[,5]/tt[,ncol(tt)] # dividing theta by length = pi
  pi=append(pi,tt)
  popage=append(popage,rep(sub(".popbams.+","",file),length(tt)))
}
pop=sub("_.+","",popage)
age=sub(".+_","",popage)

pis=data.frame(cbind(popage,pop,age,chr))
pis$pi=pi

# overall distribution of pi values
plot(density(pi),bty="n",ylab="",yaxt="n")

# pi by chromosome
ggplot(pis,aes(chr,pi))+geom_violin()


library(ggplot2)

# ggplot(pis,aes(pop,pi,fill=age))+
#   geom_violin()+
#   theme_bw()

# subtracting the per-chromosome means, for plotting

pis.centered=c()
ch="chr1"
for (ch in unique(pis$chr)){
  pchr=subset(pis,chr==ch)
  pchr$pi=pchr$pi-mean(pchr$pi)
  pis.centered=rbind(pis.centered,pchr)
}

# overall boxplot
ggplot(pis.centered,aes(pop,pi,colour=age))+
  geom_boxplot()+
  theme_bw()

# individual boxplots with per-chromosome connectors
choose="NMP"
gg=ggplot(subset(pis.centered,pop==choose),aes(age,pi))+
  geom_boxplot(outlier.shape=NA)+
  ggtitle(choose)+
  theme_bw()+
  theme(legend.position = "none")+
  geom_hline(yintercept=0,lty=3)+
  ylim(-3.1e-4,3.1e-4)+
  geom_line(aes(group=chr,color="red")) +
  geom_point(aes(group=chr),size=2,shape=21)
gg
pdf(paste(choose,"_pi.centered.pdf",sep=""),width=1.7,height=2.5)
plot(gg)
dev.off()

# deviations from mean
pis$pi.centered=pis$pi-mean(pis$pi)
ggplot(pis,aes(pop,pi.centered,colour=age))+
  geom_boxplot()
  

# computing means
library(Rmisc)
apis=data.frame(cbind(popage=rep("mean_A",14),pop=rep("mean",14),age=rep("A",14)))
apis[,4:5]=summarySE(subset(pis, age=="A"),measurevar="pi",groupvars="chr")[,c(1,3)]
jpis=data.frame(cbind(popage=rep("mean_J",14),pop=rep("mean",14),age=rep("J",14)))
jpis[,4:5]=summarySE(subset(pis, age=="J"),measurevar="pi",groupvars="chr")[,c(1,3)]
pis=rbind(pis,apis,jpis)

choose="GC"
gg=ggplot(subset(pis,pop %in% c("mean",choose)),aes(pop,pi,colour=age))+
  geom_violin()+
  ggtitle(choose)+
  theme_bw()+
  theme(legend.position = "none")

pis$pop=factor(pis$pop,levels=c("mean","GC","ST","WOR","NMP"))
ggplot(pis,aes(pop,pi,colour=age))+
  geom_boxplot()+
  
  #  geom_violin(scale="width")+
  theme_bw()+
  theme(legend.position = "none")
?geom_violin

+
  ylim(0.0038,0.0061)+
  geom_line(aes(group=chr,color="red")) +
  geom_point(aes(group=chr),size=2,shape=21)





# ----------- testing if the by-age interaction effects are significant
# linear mixed model with scalar random effect of chromosome (some chromosomes have less pi than others)

# REML method
library(lme4) 
library(lmerTest)
ll=lmer(pi~pop+pop:age+(1|chr),pis)
summary(ll)
#Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: pi ~ pop + pop:age + (1 | chr)
#    Data: pis
# 
# REML criterion at convergence: -1580.2
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.4146 -0.4267 -0.0181  0.4838  2.3395 
# 
# Random effects:
#  Groups   Name        Variance  Std.Dev. 
#  chr      (Intercept) 1.314e-07 3.625e-04
#  Residual             6.354e-09 7.971e-05
# Number of obs: 112, groups:  chr, 14
# 
# Fixed effects:
#               Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)  5.064e-03  9.920e-05  1.411e+01  51.054  < 2e-16 ***
# popNMP       5.794e-06  3.013e-05  9.100e+01   0.192  0.84792    
# popST       -8.367e-05  3.013e-05  9.100e+01  -2.777  0.00666 ** 
# popWOR      -7.345e-05  3.013e-05  9.100e+01  -2.438  0.01672 *  
# popGC:ageJ  -6.813e-05  3.013e-05  9.100e+01  -2.262  0.02611 *  
# popNMP:ageJ -1.801e-04  3.013e-05  9.100e+01  -5.979 4.34e-08 ***
# popST:ageJ   1.008e-04  3.013e-05  9.100e+01   3.347  0.00119 ** 
# popWOR:ageJ  7.834e-05  3.013e-05  9.100e+01   2.600  0.01086 *  
# ---  

#------- Bayesian method (broadly the same result as REML)
library(MCMCglmm)
mm=MCMCglmm(pi~pop+pop:age,random=~chr,data=pis,nitt=60000,burnin=10000)
HPDinterval(mm$Sol)
# lower         upper
# (Intercept)  4.852699e-03  5.271999e-03
# popNMP      -4.915871e-05  6.921060e-05
# popST       -1.448508e-04 -2.404325e-05
# popWOR      -1.343256e-04 -1.823180e-05
# popGC:ageJ  -1.279739e-04 -9.855595e-06
# popNMP:ageJ -2.432181e-04 -1.210866e-04
# popST:ageJ   3.862557e-05  1.595264e-04
# popWOR:ageJ  2.249486e-05  1.388260e-04
# attr(,"Probability")
# [1] 0.95
apply(mm$Sol,2,median)
# (Intercept)        popNMP         popST        popWOR    popGC:ageJ   popNMP:ageJ    popST:ageJ   popWOR:ageJ 
# 5.065109e-03  5.753827e-06 -8.395416e-05 -7.411881e-05 -6.874904e-05 -1.810929e-04  1.007064e-04  7.963282e-05 


#-------- SFS of adults vs juveniles

infiles=c("GC.sfs","nmp.sfs","st.sfs","wor.sfs")
source("sfs2matrix.R")
filename=infiles[2]

sfs2matrix=function(sfs,n1,n2) {
  dd=matrix(ncol=2*n1+1,nrow=2*n2+1,sfs)
    dd[1,1]=dd[2*n2+1,2*n1+1]=0
  return(apply(dd,2,rev))
}

plotSFS=function(filename,n1,n2){
  require(RColorBrewer)
  require(pheatmap)
  s1d=t(read.table(filename))[,1]
  ss=log(sfs2matrix(s1d,n1,n2)+1e-2,10)
  col=colorRampPalette(brewer.pal(n = 7, name = "Greys"))(100)
  pp=pheatmap(ss,cluster_rows = F, cluster_cols = F,border_color = F,color=col,angle_col=0)
  return(pp)
}

p=plotSFS("GC.sfs",33,20)
p=plotSFS("nmp.sfs",31,34)
p=plotSFS("wor.sfs",30,37)
p=plotSFS("st.sfs",32,43)
