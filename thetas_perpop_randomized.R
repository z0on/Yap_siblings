setwd("~/Dropbox/yap_hetcheck_2020/Yap_siblings") 
system("ls *popbamsR.pestPG >ins")
infiles=scan("ins", what="c")
library(ggplot2)

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
choose="WOR"
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
pdf(paste(choose,"_pi.centeredR.pdf",sep=""),width=1.7,height=2.5)
plot(gg)
dev.off()


# ----------- testing if the by-age interaction effects are significant
# linear mixed model with scalar random effect of chromosome (some chromosomes have less pi than others)

# REML method
library(lme4) 
library(lmerTest)
ll=lmer(pi~pop+pop:age+(1|chr),pis)
summary(ll)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: pi ~ pop + pop:age + (1 | chr)
#    Data: pis
# 
# REML criterion at convergence: -1580
# 
# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -2.47216 -0.56062 -0.00493  0.58671  2.35672 
# 
# Random effects:
#  Groups   Name        Variance  Std.Dev. 
#  chr      (Intercept) 1.306e-07 3.614e-04
#  Residual             6.369e-09 7.981e-05
# Number of obs: 112, groups:  chr, 14
# 
# Fixed effects:
#               Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)  5.016e-03  9.892e-05  1.412e+01  50.706   <2e-16 ***
# popNMP       2.995e-06  3.016e-05  9.100e+01   0.099   0.9211    
# popST       -9.790e-06  3.016e-05  9.100e+01  -0.325   0.7463    
# popWOR       2.429e-05  3.016e-05  9.100e+01   0.805   0.4228    
# popGC:ageJ  -1.303e-05  3.016e-05  9.100e+01  -0.432   0.6669    
# popNMP:ageJ -1.383e-05  3.016e-05  9.100e+01  -0.458   0.6477    
# popST:ageJ   5.502e-05  3.016e-05  9.100e+01   1.824   0.0714 .  
# popWOR:ageJ -1.766e-05  3.016e-05  9.100e+01  -0.585   0.5597    

#------- Bayesian method (broadly the same result as REML)
library(MCMCglmm)
mm=MCMCglmm(pi~pop+pop:age,random=~chr,data=pis,nitt=60000,burnin=10000)
HPDinterval(mm$Sol)
# lower        upper
# (Intercept)  4.795829e-03 5.220395e-03
# popNMP      -5.170403e-05 6.751507e-05
# popST       -6.813658e-05 5.155490e-05
# popWOR      -3.487517e-05 8.294893e-05
# popGC:ageJ  -6.982390e-05 4.960284e-05
# popNMP:ageJ -7.484647e-05 4.563459e-05
# popST:ageJ  -6.103894e-06 1.149283e-04
# popWOR:ageJ -7.822961e-05 3.879889e-05
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
