setwd("~/Dropbox/yap_hetcheck_2020/Yap_siblings") # change this to where your scp'd files are

#------ IBS of all bams

allbams=scan("bams.qc",what="c")
allbams=sub(".bam","",allbams)
ibsma=as.matrix(read.table("g3all.ibsMat"))
dimnames(ibsma)=list(allbams,allbams)
hc=hclust(as.dist(ibsma),method="ave")
plot(hc,cex=0.0001)


#---- removing one per sib

togo=c("NMP_143J","NMP_154J")

#------- removing all but one clone representative (keeping the best-covered one)

plot(hc,cex=0.5)
abline(h=0.125,col="red")
# looks like if we cut the tree here (at 0.05), groups of 2-4 individuals will be clonal groups
cuts=cutree(hc,h=0.125)
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

goods=allbams[allbams %in% c(uniqs,bestclones)]
goods=goods[!(goods %in% togo)]
length(goods)

ibsma1=ibsma[goods,goods]
hc=hclust(as.dist(ibsma1),method="ave")
plot(hc,cex=0.5)

cc=capscale(as.dist(ibsma1)~1)
plot(cc)
hetouts=c("GC_29J","WOR_10A")
scc=data.frame(scores(cc,scaling=1)$sites)
points(scc[row.names(scc) %in% hetouts,],col="red",pch=19)
#----- writing good bams list

length(goods)
# 261
goods=paste(goods,".bam",sep="")
write(goods,file="bams2")
