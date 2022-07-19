setwd("~/Dropbox/yap_hetcheck_2020/") # change this to where your scp'd files are
library(pheatmap)

# reading long relatedness table, output of ngsRelate
rel=read.table("g3.relatedness",sep="\t",header=T)

# creating an empty square matrix
relm=matrix(0,nrow=length(unique(rel$a))+1,ncol=length(unique(rel$a))+1)

library(tidyverse)
# filling up the square matrix with entries from "rab" column
sfs=list()
for (a in unique(c(rel$a))) {
	for (b in unique(rel$b)) {
		if (b<=a) {next}
	  relm[a+1,b+1]=relm[b+1,a+1]=rel[rel$a==a & rel$b==b,"rab"]
	}
}
diag(relm)=1
dim(relm)

relm[sibs,sibs]

# adding names to columns and rows - assuming ngsRelate was run on angsd result obtained for bams.nr file
bams=scan("goodbams.sub100k",what="character") # list of bam files
bams=sub(".sub100k.bam","",bams)
dimnames(relm)=list(bams,bams)
bams[grep("NMP_147",bams)]

# reading sfs
sfs=list();i=1;sfsnames=c()
for (a in unique(c(rel$a))) {
  for (b in unique(rel$b)) {
    if (b<=a) {next}
    mm=as.numeric(strsplit(rel[rel$a==a & rel$b==b,"X2dsfs"],",")[[1]])
#    sfs[[i]]=as.matrix(rbind(mm[7:9],mm[4:6],mm[1:3]))
    sfs[[i]]=mm
        i=i+1
    sfsnames=c(sfsnames,paste(bams[a+1],bams[b+1],sep="_"))
  }
}


# reading inbreedings
inb=c()
for (a in unique(c(rel$a))) {
    inb=append(inb,mean(as.numeric(rel[rel$a==a,"Fa"])))
}
inb=append(inb,mean(as.numeric(rel[rel$b=="261","Fa"])))
names(inb)=bams

source("~/Dropbox/sfs2matrix.R")
plotSFS=function(pairname){
  s1d=sfs[[which(sfsnames==pairname)]]
  ss=log(sfs2matrix(s1d,1,1)+1e-2,10)
  dimnames(ss)=list(rev(c("AA","Aa","aa")),c("AA","Aa","aa"))
  col=colorRampPalette(brewer.pal(n = 7, name = "Greys"))(100)
  pp=pheatmap(ss,cluster_rows = F, cluster_cols = F,border_color = F,color=col,angle_col=0)
  plot(pp)
}
10^(-1.5)
10^(-0.5)

pdf("NMP_139J_NMP_143J_sfs.pdf",height=1.5,width=2)
plotSFS("NMP_139J_NMP_143J")
dev.off()

pdf("NMP_147J_NMP_154J_sfs.pdf",height=1.5,width=2)
plotSFS("NMP_147J_NMP_154J")
dev.off()

pdf("randomPair_sfs.pdf",height=1.5,width=2)
plotSFS(sample(sfsnames,1))
dev.off()

# plotting the matrix
library(pheatmap)
pheatmap(relm,cex=0.5)


# removing pair of sibs to do capscale
sibs2go=which(bams %in% c("NMP_154J","NMP_143J"))
relm.nosibs=relm[-sibs2go,-sibs2go]
reld=as.dist(1-relm.nosibs)
library(vegan)
rcap=capscale(reld~1)

#predicting full dataset (with sibs)
pred=predict(rcap,relm,type='sp',scaling="sites")
rscor=data.frame(scores(pred,scaling=1)[,1:2])

rscor$inbreed=inb
sibnames=c("NMP_147J","NMP_154J","NMP_139J","NMP_143J")

i2p=data.frame(row.names(rscor),sub("_\\d+","_",row.names(rscor)))
i2p[,2]=sub("_H.+","_A",i2p[,2])
write.table(i2p,quote=F, col.names=F,row.names=F,file="inds2pops")

# tracking nimpal juveniles
nmpj=rep("other",nrow(rscor))
nmpj[grep("NMP_\\d+J",row.names(rscor))]="NMP Juv"

# tracking sibs
sibs1=sibs2=rep(FALSE,nrow(rscor))
sibs1[which(row.names(rscor) %in% c("NMP_147J","NMP_154J"))]=TRUE
sibs2[which(row.names(rscor) %in% c("NMP_139J","NMP_143J"))]=TRUE
rscor$nmp=nmpj

rscor$pcinbreed=read.table("pcangsd.min.inbreed.samples")[,1]
plot(pcinbreed~inbreed,rscor)
load("pcangsd_admixture.RData")
rscor$adm=tbl$V2

pdf("pcoa_relatedness_inbreedScale.pdf",height=4,width=4)
ggplot()+geom_point(rscor,mapping=aes(MDS1,MDS2,color=nmpj,size=inbreed),pch=1)+theme_bw()+coord_equal()+
  geom_point(rscor[sibs1,],mapping=aes(MDS1,MDS2),pch=3,color="black")+
  geom_point(rscor[sibs2,],mapping=aes(MDS1,MDS2),pch=4,color="black")
dev.off()

# loading admixture, setting admixture values for missing sibs to the same as the other sib
load("pcangsd_admixture.RData")
head(tbl)
rscor$admix=tbl[row.names(rscor),"V2"]

pdf("pcoa_relatedness_admixture.pdf",height=4,width=4)
ggplot()+geom_point(rscor,mapping=aes(MDS1,MDS2,color=nmpj,size=admix),pch=1)+
  scale_size_continuous(range=c(0.5,5),trans="exp")+
  theme_bw()+coord_equal()+
  geom_point(rscor[sibs1,],mapping=aes(MDS1,MDS2),pch=3,color="black")+
  geom_point(rscor[sibs2,],mapping=aes(MDS1,MDS2),pch=4,color="black")
dev.off()

#correlation with admixture (perfect one)
summary(lm(MDS1~admix,rscor[!(row.names(rscor) %in% c("NMP_143J","NMP_154J")),]))
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.135990   0.008627   15.76   <2e-16 ***
#   admix       -0.272962   0.016195  -16.86   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05333 on 258 degrees of freedom
# Multiple R-squared:  0.5241,	Adjusted R-squared:  0.5222
plot(MDS1~admix,rscor)

rscor$pop_age=factor(tbl$pop)
rscor$pop=factor(sub("_[AJ]","",rscor$pop_age))
rscor$age=factor(sub(".+_","",rscor$pop_age))

# heterozygosity
zy=read.table("goodbams.minsites.het")
names(zy)=c("sample","nsites","hz")
zy$sample=sub("\\..+","",zy$sample)
row.names(zy)=zy$sample
zy$lns=log(zy$nsites,10)
plot(hz~lns,zy)
abline(lm(hz~lns,zy))

# bams=scan("bams.sub100k","c")
# bams=sub("\\..+","",bams)
# zy=data.frame(cbind(sample=bams,hz=SFS[2,]))
# row.names(zy)=bams
rscor$hz=as.numeric(zy[row.names(rscor),"hz"])
rscor$lns=as.numeric(zy[row.names(rscor),"lns"])

head(rscor)
save(rscor,file="relatedness_pcoa_admix_inbreed_het_nsites_MinSites.RData")
rscor[sibs,]
pairs(rscor[,c(1:3,7,11,12)])

summary(lm(rscor[,1]~nmpj))
# NS

summary(lm(MDS1~inbreed,rscor))
# NS
summary(lm(rscor$inbreed~nmpj))
# p=0.015

pdf("pcoa_relatedness_heterozygosity.pdf",height=2.5,width=3.8)
ggplot()+geom_point(rscor,mapping=aes(MDS1,MDS2,color=nmpj,size=hz),pch=1)+
#  scale_size_continuous(range=c(0.5,5),trans="exp")+
  theme_bw()+coord_equal()+
  geom_point(rscor[sibs1,],mapping=aes(MDS1,MDS2),pch=3,color="black")+
  geom_point(rscor[sibs2,],mapping=aes(MDS1,MDS2),pch=4,color="black")
dev.off()

ggplot()+geom_point(rscor,mapping=aes(MDS1,MDS2,color=nmpj,size=lns),pch=1)+
  #  scale_size_continuous(range=c(0.5,5),trans="exp")+
  theme_bw()+coord_equal()+
  geom_point(rscor[sibs1,],mapping=aes(MDS1,MDS2),pch=3,color="black")+
  geom_point(rscor[sibs2,],mapping=aes(MDS1,MDS2),pch=4,color="black")

plot(MDS1~lns,rscor)
plot(admix~lns,rscor)

plot(inbreed~hz,rscor)
rscor[sibs,]
plot(rscor$hz[order(rscor$hz)])
plot(rscor$lns[order(rscor$lns)])
rscor[order(rscor$hz),]
rscor[order(rscor$lns),]


#------- IBS ordination

ibsm=as.matrix(read.table("g3.ibsMat"))
dimnames(ibsm)=list(bams,bams)

ibsm.nosibs=relm[-sibs2go,-sibs2go]
ibsd=as.dist(1-ibsm.nosibs)
library(vegan)
ibscap=capscale(ibsd~1)

#predicting full dataset (with sibs)
predi=predict(ibscap,1-ibsm,type='sp',scaling="sites")
ibsscor=rscor
ibsscor[,1:2]=scores(predi,scaling=1)[,1:2]
save(ibsscor,file="IBSscores_etc.RData")


#-------plotting ordination

plotOrdi=function(ordination,groupingFactor,spider=T,ellipse=T,...){
  require(adegenet) # for transp()
  palette(rainbow(length(levels(groupingFactor))))
  colors=as.numeric(groupingFactor)
  colpops=as.numeric(as.factor(sort(levels(groupingFactor))))
  plot(ordination,type="n",mgp=c(2.3,1,0))
  points(ordination,pch=19,col=transp(colors,alpha=0.5))
  if(spider==T) { ordispider(ordination,groupingFactor,label=T,col=colpops,... )}
  if(ellipse==T) { ordiellipse(ordination,groupingFactor,label=T,col=colpops, draw="polygon",...)}
}

par(mfrow=c(1,1))
pdf("ordination_pop_age.pdf",height=4.8,width=4)
plotOrdi(ibsscor[,1:2],rscor$pop_age,ellipse=F,cex=0.8)
dev.off()

pdf("ordination_pop.pdf",height=4.8,width=4)
plotOrdi(predi,rscor$pop,ellipse=F,cex=0.8)
dev.off()

pdf("ordination_age.pdf",height=4.8,width=4)
plotOrdi(pred,rscor$age,ellipse=F,cex=0.8)
dev.off()

#-------------------- plotting IBS with sibs

pdf("pcoa_IBS_ellipse_sibs.pdf",height=3,width=3.5)
ggplot()+geom_point(ibsscor,mapping=aes(MDS1,MDS2,color=pop_age),pch=1)+
  #  scale_size_continuous(range=c(0.5,5),trans="exp")+
  theme_void()+coord_equal()+stat_ellipse(ibsscor,mapping=aes(MDS1,MDS2,color=pop_age),type="norm")+
  geom_point(ibsscor[sibs1,],mapping=aes(MDS1,MDS2),pch=3,size=5,color="grey20")+
  geom_point(ibsscor[sibs2,],mapping=aes(MDS1,MDS2),pch=4,size=5,color="grey20")
dev.off()

ibsscor[sibs,]
boxplot(ibsscor$hz)
abline(h=0.00549239)

# NMP_143J is a heterozygosity outlier, but there are two more like him and worse

# by pop

adonis(pred~pop+age,data=rscor)
adonis(pred~pop_age,data=rscor)

ordiellipse(pred,rscor$pop,label=T)
ordispider(pred,rscor$pop,label=T)
adonis()


