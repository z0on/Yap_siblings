setwd("~/Dropbox/yap_hetcheck_2020/Yap_siblings") # change this to where your scp'd files are

#------ IBS of all bams

allbams=scan("goodbams.nohh",what="c")
allbams=sub(".bam","",allbams)
ibsma=as.matrix(read.table("all.nohh.ibsMat"))
dimnames(ibsma)=list(allbams,allbams)
hc=hclust(as.dist(ibsma),method="ave")
pdf("allbams_hclust.pdf",height=3.5,width=8)
plot(hc,cex=0.0001)
dev.off()

#---- removing wrong collections

plot(hc,cex=0.5)
abline(h=0.1,col="red")
# looks like if we cut the tree here (at 0.1), the larger group will be our species and the smaller one - wrong species

cuts=cutree(hc,h=0.1)
goodsp=names(cuts)[cuts==1]

#------- removing all but one clone representative (keeping the best-covered one)

plot(hc,cex=0.5)
abline(h=0.05,col="red")
# looks like if we cut the tree here (at 0.05), groups of 2-4 individuals will be clonal groups
cuts=cutree(hc,h=0.05)
cutgroups=table(cuts)
clones=cutgroups[cutgroups>1]
uni=cutgroups[cutgroups==1]
uniqs=names(cuts[cuts %in% names(uni)])

# reading qualities
qu=read.table("quality.txt")
row.names(qu)=gsub(".trim.bt2.bam|/work/03121/sbarfie/yap/|.bam","",qu$V1)

# selecting best clone representatives
bestclones=c()
for (cl in names(clones)){
  cls=names(cuts[cuts==cl])
   clqu=qu[cls,]$V2
  clrep=cls[which(clqu==max(clqu)[1])]
  bestclones=c(bestclones,clrep)
}

goods=goodsp[goodsp %in% c(uniqs,bestclones)]

#----- sanity check: plotting hctree after removal of bads

ibsma1=ibsma[goods,goods]
hc=hclust(as.dist(ibsma1),method="ave")
plot(hc,cex=0.5)

#----- writing good bams list

length(goods)
# 261
goods=paste(goods,".bam",sep="")
write(goods,file="bams2")
