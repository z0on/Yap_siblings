source('~/Dropbox/yap_hetcheck_2020/plot_admixture_v5_function.R')
library(RcppCNPy)

# assembling the input table
dir="~/Dropbox/yap_hetcheck_2020/Yap_siblings/"
setwd(dir)
tbl=read.table("pcangsd.min.admix.2.Q")
npops=ncol(tbl)

pops=sub("\\d+","",bams)
pops=sub("_H","_A",pops)
i2p=data.frame(cbind(bams,pops))
names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
row.names(tbl)=tbl$ind
tbl$pop=factor(tbl$pop)

save(tbl,file="pcangsd_admixture.RData")

# saving 8 least-admixed blue group guys and 8 least-admixed blue group guys for SFS and theta analysis
tops=row.names(tbl[order(tbl$V2,decreasing=T),])[1:8] # all NMP juveniles
bottoms=row.names(tbl[order(tbl$V2,decreasing=F),])[2:9] # all WOR adults
write(paste(tops,".sub100k.bam",sep=""),file="bluetops.bams")
write(paste(bottoms,".sub100k.bam",sep=""),file="redtops.bams")

pdf("pcangsd_admix.pdf",height=2.5,width=8)
ords=plotAdmixture(data=tbl,npops=npops,vshift=0.05,angle=45,cex=0.8,grouping.method="sequential")
dev.off()

