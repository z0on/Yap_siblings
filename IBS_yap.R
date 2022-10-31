setwd("~/Dropbox/yap_hetcheck_2020/Yap_siblings") # change this to where your scp'd files are
library(pheatmap)

#------ IBS of all bams

allbams=scan("goodbams.nohh",what="c")
allbams=sub(".bam","",allbams)
ibsma=as.matrix(read.table("all.nohh.ibsMat"))
dimnames(ibsma)=list(allbams,allbams)
hc=hclust(as.dist(ibsma),method="ave")
pdf("allbams_hclust.pdf",height=4,width=8)
plot(hc,cex=0.0001)
dev.off()
plot(hc,cex=0.5)

# ---- list of bams without clones, wrong collections, and heterozygosity outliers

bams=scan("goodbams.sub100k",what="c")
bams=sub(".sub100k.bam","",bams)

#---- quality (coverage) of sibs

qq=read.table("quality.txt",header=F)
names(qq)=c("sample","quality")
qq$sample=sub("/work/03121/sbarfie/yap/","",qq$sample)
qq$sample=sub(".trim.bt2.bam","",qq$sample)
row.names(qq)=qq$sample
qq=qq[bams,]
qq=qq[order(qq$quality),]
sibs=c("NMP_139J", "NMP_143J", "NMP_147J", "NMP_154J")
plot(qq$quality,xaxt="n",bty="n", xlab="",ylab="quality",mgp=c(2.3,1,0),type="l")
#boxplot(qq$quality,ylab="quality")
sibqq=qq[sibs,"quality"]
abline(h=sibqq,lty=3,col="red")

#---- nsites of sibs

sh=read.table("goodbams.minsites.het")
names(sh)=c("sample","nsites","het")
sh$sample=sub(".bam","",sh$sample)
qqsh=merge(qq,sh,by="sample")
head(qqsh)
plot(nsites~quality,qqsh,log="y")

pdf("sibQuality.pdf",height=4,width=3.3)
plot(nsites~quality,qqsh,col="grey70",mgp=c(2.3,1,0))
points(nsites~quality,qqsh[qqsh$sample %in% c("NMP_147J", "NMP_154J"),],pch=3,cex=2)
points(nsites~quality,qqsh[qqsh$sample %in% c("NMP_139J", "NMP_143J"),],pch=4,cex=2)
#points(nsites~quality,qqsh[qqsh$sample %in% c("NMP_143J"),],pch=16)
dev.off()

#------- IBS ordination
ibsm=as.matrix(read.table("g3.ibsMat"))
dimnames(ibsm)=list(bams,bams)

# removing pair of sibs to do capscale
sibs2go=which(bams %in% c("NMP_154J","NMP_143J"))
ibsm.nosibs=ibsm[-sibs2go,-sibs2go]
ibsd=as.dist(1-ibsm.nosibs)
library(vegan)
ibscap=capscale(ibsd~1)

#predicting full dataset (with sibs)
predi=predict(ibscap,1-ibsm,type='sp',scaling="sites")
ibsscor=data.frame(scores(predi,scaling=1)[,1:2])

sibnames=c("NMP_147J","NMP_154J","NMP_139J","NMP_143J")

# tracking nimpal juveniles
nmpj=rep("other",nrow(ibsscor))
nmpj[grep("NMP_\\d+J",row.names(ibsscor))]="NMP Juv"

# tracking sibs
sibs1=sibs2=rep(FALSE,nrow(ibsscor))
sibs1[which(row.names(ibsscor) %in% c("NMP_147J","NMP_154J"))]=TRUE
sibs2[which(row.names(ibsscor) %in% c("NMP_139J","NMP_143J"))]=TRUE
ibsscor$nmp=nmpj

ibsscor$pcinbreed=read.table("pcangsd.min.inbreed.samples")[,1]

# loading admixture, setting admixture values for missing sibs to the same as the other sib
load("pcangsd_admixture.RData")
head(tbl)
ibsscor$admix=tbl[row.names(ibsscor),"V2"]

ibsscor$pop_age=factor(tbl$pop)
ibsscor$pop=factor(sub("_[AJ]","",ibsscor$pop_age))
ibsscor$age=factor(sub(".+_","",ibsscor$pop_age))

# heterozygosity

zy=read.table("goodbams.minsites.het")
names(zy)=c("sample","nsites","hz")
zy$sample=sub("\\..+","",zy$sample)
row.names(zy)=zy$sample
zy$lns=log(zy$nsites,10)
plot(hz~lns,zy)
abline(lm(hz~lns,zy))

ibsscor$hz=as.numeric(zy[row.names(ibsscor),"hz"])
ibsscor$lns=as.numeric(zy[row.names(ibsscor),"lns"])

head(ibsscor)
save(ibsscor,file="IBS_pcoa_admix_pcinbreed_het_nsites_MinSites.RData")

pairs(ibsscor[,c(1:3,7,11,12)])

summary(lm(ibsscor[,1]~nmpj))
# NS

summary(lm(MDS1~inbreed,ibsscor))
# NS
summary(lm(ibsscor$inbreed~nmpj))
# p=0.015

#----- heterozygosity/nsites outliers removal

ll=load("IBS_pcoa_admix_pcinbreed_het_nsites_MinSites.RData")
dim(ibsscor)
plot(density(ibsscor$hz))
plot(density(ibsscor$lns))
ibsscor$hz_z=abs(scale(ibsscor$hz))
ibsscor$lns_z=abs(scale(ibsscor$lns))
table(ibsscor$hz_z>2)
# FALSE  TRUE 
# 251    11
table(ibsscor$lns_z>2)
# FALSE  TRUE 
# 256     6
table(ibsscor$lns_z>2 | ibsscor$hz_z>2)
# FALSE  TRUE 
# 250    12 

pass=ibsscor[(ibsscor$lns_z<2 & ibsscor$hz_z<2),]

# writin down per-pop per-cohort bam lists
for(pa in levels(ibsscor$pop_age)) {
  bb=row.names(ibsscor[ibsscor$pop_age == pa,])
  write(paste(bb,".bam",sep=""),file=paste(pa,".bams",sep=""))
  print(paste(pa,length(bb)))
}
# [1] "GC_A 33"
# [1] "GC_J 21"
# [1] "NMP_A 31"
# [1] "NMP_J 34"
# [1] "ST_A 32"
# [1] "ST_J 43"
# [1] "WOR_A 31"
# [1] "WOR_J 37"

#-------------------- plotting IBS with sibs

pdf("pcoa_IBS_ellipse075_sibs.pdf",height=3,width=3.5)
ggplot()+geom_point(ibsscor,mapping=aes(MDS1,MDS2,color=pop_age),pch=1)+
  #  scale_size_continuous(range=c(0.5,5),trans="exp")+
  theme_void()+coord_equal()+stat_ellipse(ibsscor,mapping=aes(MDS1,MDS2,color=pop_age),type="norm",level=0.75)+
  geom_point(ibsscor[sibs1,],mapping=aes(MDS1,MDS2),pch=3,size=5,color="grey20")+
  geom_point(ibsscor[sibs2,],mapping=aes(MDS1,MDS2),pch=4,size=5,color="grey20")
dev.off()

# detecting structure
adonis2(rscor[,1:2]~pop_age,data=rscor,method="euclidean")
# Df SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)   
# pop_age     7 0.0009001 1.2859e-04   2.085 0.05434  0.006 **
#   Residuals 254 0.0156647 6.1672e-05         0.94566          
# Total     261 0.0165648                    1.00000   


# is it because of WOR_A? YES

POP="WOR_A";ps=c()
for(POP in levels(ibsscor$pop_age)) {
  selpop=rep(FALSE,nrow(ibsscor))
  selpop[ibsscor$pop_age==POP]=TRUE
  selpop=data.frame(selpop)
  aa=adonis(ibsscor[,1:2]~selpop,data=selpop,method="eucl",permutations=999)
  ps=c(ps,aa$aov.tab$`Pr(>F)`[1])
}
data.frame(cbind(pop_age=levels(ibsscor$pop_age),pval=ps))

adonis(subset(ibsscor,pop_age != "WOR_A")[,1:2]~pop_age,data=subset(ibsscor,pop_age != "WOR_A"),method="euclidean")
# Df SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)
# pop_age     6 0.0002683 4.4723e-05 0.75359 0.01979  0.692
# Residuals 224 0.0132938 5.9347e-05         0.98021       
# Total     230 0.0135621                    1.00000  



?stat_ellipse

#-------plotting ordination vegan way

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
plotOrdi(rscor[,1:2],rscor$pop_age,ellipse=F,cex=0.8)
dev.off()

pdf("ordination_pop.pdf",height=4.8,width=4)
plotOrdi(rscor,rscor$pop,ellipse=F,cex=0.8)
dev.off()

pdf("ordination_age.pdf",height=4.8,width=4)
plotOrdi(pred,ibsscor$age,ellipse=F,cex=0.8)
dev.off()

