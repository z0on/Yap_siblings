setwd("~/Dropbox/yap_hetcheck_2020/Yap_siblings") # change this to where your scp'd files are
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

# adding names to columns and rows - assuming ngsRelate was run on angsd result obtained for bams2 file
bams=scan("bams2",what="character") # list of bam files
bams=sub(".bam","",bams)
dimnames(relm)=list(bams,bams)

# reading pairwise sfs (from ngsRelate ourput)
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

# reading inbreeding (from ngsRelate ourput)

rel$aa=paste("x",rel$a)
ggplot(rel,aes(aa,Fa))+geom_boxplot()
rel=na.omit(rel)
inb=c()
for (a in unique(rel$a)) {
    inb=append(inb,mean(as.numeric(rel[rel$a==a,"Fa"])))
}
inb=append(inb,mean(as.numeric(rel[rel$b=="260","Fb"])))
names(inb)=bams
table(is.na(inb))

# plotting AFS between sibs and between random samples

source("sfs2matrix.R")
plotSFS=function(pairname){
  require(RColorBrewer)
  require(pheatmap)
  s1d=sfs[[which(sfsnames==pairname)]]
  ss=log(sfs2matrix(s1d,1,1)+1e-2,10)
  dimnames(ss)=list(rev(c("AA","Aa","aa")),c("AA","Aa","aa"))
  col=colorRampPalette(brewer.pal(n = 7, name = "Greys"))(100)
  pp=pheatmap(ss,cluster_rows = F, cluster_cols = F,border_color = F,color=col,angle_col=0)
  plot(pp)
}

# pdf("NMP_139J_NMP_143J_sfs.pdf",height=1.5,width=2)
# plotSFS("NMP_139J_NMP_143J")
# dev.off()

pdf("NMP_147J_NMP_154J_sfs.pdf",height=1.5,width=2)
plotSFS("NMP_147J_NMP_154J")
dev.off()

pdf("randomPair_sfs.pdf",height=1.5,width=2)
plotSFS(sample(sfsnames,1))
dev.off()

#----- plotting hclust tree based on relatedness

rhc=hclust(as.dist(1-relm),method="ave")
pdf("relatedness_hclust.pdf",height=3.5,width=8)
plot(rhc,cex=0.001)
dev.off()

#------- IBS ordination

ibsm=as.matrix(read.table("g3.ibsMat"))
#ibsm=cov2cor(as.matrix(read.table("pcangsd.min.cov")))
dimnames(ibsm)=list(bams,bams)

# removing one of the sibs to build ordination (we will put him back later)
sibs2go="NMP_154J"
sibs1=c("NMP_154J","NMP_147J")
ibsm.nosibs=ibsm[bams!="NMP_154J",bams!="NMP_154J"]
ibsd=as.dist(1-ibsm.nosibs)
library(vegan)
ibscap=capscale(ibsd~1)
plot(ibscap$CA$eig/sum(ibscap$CA$eig))

#predicting full dataset (with the sib that has been withheld)
predi=predict(ibscap,1-ibsm,type='sp',scaling="sites")
ibsscor=data.frame(scores(predi,scaling=1)[,1:2])
ibsscor$inbreed=inb
ibsscor$pop=sub("_.+","",row.names(ibsscor))
ibsscor$age="A"
ibsscor$age[grep("J",row.names(ibsscor))]="J"
ibsscor$pop_age=paste(ibsscor$pop,ibsscor$age,sep="_")

pdf("pcoa_IBS_ellipse_sibs.pdf",height=3,width=3.5)
ggplot()+geom_point(ibsscor,mapping=aes(MDS1,MDS2,color=pop_age),pch=1)+
  #  scale_size_continuous(range=c(0.5,5),trans="exp")+
  theme_void()+coord_equal()+stat_ellipse(ibsscor,mapping=aes(MDS1,MDS2,color=pop_age),type="norm",level=0.75)+
  geom_point(ibsscor[sibs1,],mapping=aes(MDS1,MDS2),pch=3,size=5,color="grey20")
#  geom_point(ibsscor[sibs2,],mapping=aes(MDS1,MDS2),pch=4,size=5,color="grey20")
dev.off()

#------ reading quailty
qq=read.table("quality.txt",sep=" ")
qq$V1=sub(".+yap.","",qq$V1)
qq$V1=sub("\\..+","",qq$V1)
row.names(qq)=qq$V1
qq$V1=NULL
ibsscor$qual=qq[row.names(ibsscor),"V2"]

# ------ reading heterozygosity and nsites
zz=read.table("goodsites.het")
row.names(zz)=sub("\\.bam","",zz$V1)
ibsscor$log.nsites=log(zz[row.names(ibsscor),"V2"],10)
ibsscor$het=zz[row.names(ibsscor),"V3"]

#-------- reading pcangsd assignments and inbreeding

pa=read.table("pcangsd.min.admix.2.Q")
row.names(pa)=bams
ibsscor$q1=pa[row.names(ibsscor),"V1"]
ibsscor$q2=pa[row.names(ibsscor),"V2"]

pai=read.table("pcangsd.min.inbreed.samples")
ibsscor$pa.inb=pai$V1

pairs(ibsscor[,c("qual","inbreed","pa.inb","log.nsites","het","q1","MDS1")],cex=0.5,col=rgb(0,0,0,alpha=0.5))

# three "site-standard" samples show outlyingly high inbreeding - artifact! ignore them in plotting below
ignor=row.names(head(ibsscor[order(ibsscor$pa.inb,decreasing=T),],3))
ibsscor=ibsscor[!(row.names(ibsscor) %in% ignor),]
dim(ibsscor)

# ------ plotting inbreeding variation
ibsscor$pop=factor(ibsscor$pop,levels=c("GC","ST","WOR","NMP"))
ibsscor$pa.inb.c=ibsscor$pa.inb-mean(ibsscor$pa.inb)
gg=ggplot(ibsscor,aes(pop,pa.inb.c,color=age))+
  geom_boxplot()+theme_bw()+
  geom_hline(yintercept=0,lty=3)+
  theme(legend.position = "bottom")+xlab("")+ylab("inbr.dev.")


pdf("inbreed.pdf",height=2.5,width=2.5)
plot(gg)
dev.off()
#------ regressing quality out of inbreeding

plot(inbreed~qual,ibsscor)
ibsscor=na.omit(ibsscor)
library(MASS) # for rlm, robust regression
ibsscor$inb=residuals(rlm(inbreed~qual,ibsscor))
plot(inb~qual,ibsscor)

plot(inbreed~MDS1,ibsscor)
plot(inbreed~MDS2,ibsscor)

head(ibsscor[order(ibsscor$inbreed,decreasing=T),])

ggplot(ibsscor,aes(pop,inb,color=age))+geom_boxplot()

# -----detecting structure in IBS data

adonis2(as.dist(1-ibsm)~pop_age,data=ibsscor,method="euclidean",permutations=9999)
screeplot(ibscap)
adonis2(ibsscor[,1:2]~pop_age,data=ibsscor,method="euclidean",permutations=9999)
# # Df  SumOfSqs      R2    F Pr(>F)  
# Df  SumOfSqs      R2      F Pr(>F)
# pop_age    7 0.0001381 0.02626 0.9747 0.4848
# Residual 253 0.0051195 0.97374              
# Total    260 0.0052576 1.00000   


# checking whether each pop_age group is significantly different from the rest
POP="WOR_A";ps=c()
for(POP in unique(ibsscor$pop_age)) {
  selpop=rep(FALSE,nrow(ibsscor))
  selpop[ibsscor$pop_age==POP]=TRUE
  selpop=data.frame(selpop)
  aa=adonis2(ibsscor[,1:2]~selpop,data=selpop,method="eucl",permutations=999)
  ps=c(ps,aa$`Pr(>F)`[1])
}
data.frame(cbind(pop_age=unique(ibsscor$pop_age),pval=ps))
# pop_age  pval
# 1    GC_A 0.514
# 2    GC_J 0.159
# 3   NMP_A 0.731
# 4   NMP_J  0.98
# 5    ST_J 0.535
# 6    ST_A 0.118
# 7   WOR_A 0.496
# 8   WOR_J 0.212
# 

#------ writing bam lists for all pop_age groups

for (p in unique(ibsscor$pop_age)) {
 bb=row.names(subset(ibsscor,pop_age==p))
 fname=paste(p,".popbams",sep="")
 write(paste(bb,".bam",sep=""),fname)
}  

# writing randomized bam groups
ibsscor$pop_age_random=sample(ibsscor$pop_age)
for (p in unique(ibsscor$pop_age)) {
  bb=row.names(subset(ibsscor,pop_age_random==p))
  fname=paste(p,".popbamsR",sep="")
  write(paste(bb,".bam",sep=""),fname)
}  

#----------- plottign supplemental figures
sibs1=c("NMP_154J","NMP_147J")
plot(log.nsites~qual,ibsscor,xlab="depth-quality",ylab="log10(N sites)",col=rgb(0,0,0,alpha=0.3),cex=0.8,mgp=c(2.3,1,0))
points(log.nsites~qual,ibsscor[sibs1,],pch=3,col="red",cex=1.6)

plot(het~qual,ibsscor,xlab="depth-quality",ylab="heterozygosity",col=rgb(0,0,0,alpha=0.3),cex=0.8,mgp=c(2.3,1,0))
points(het~qual,ibsscor[sibs1,],pch=3,col="red",cex=1.6)

