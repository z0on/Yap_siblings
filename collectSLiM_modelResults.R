
# repace directory location below with the name of the cloned github repo directory
setwd("~/Dropbox/yap_hetcheck_2020")

modelRuns=scan("modnames",what="c")
modelRuns=c("norm")
load('sst_modified.RData')

library(tidyverse)
library(Rmisc)
run=c("a","b","c","d")
for (s in modelRuns ){
	respranges=list();prewarms=list();raws=list();stds=list();relresps=list();hs=list()
	r="a";sets=s;index=""
	message(sets)
	for (r in run){
		message(paste("   ",r))
	  system(paste("cat", paste(sets,r,sep=""),"| perl -pe 's/^no/#no/' | perl -pe 's/^empty/#empty/'>temp"))
	   g85=read.table("temp",sep="\t",skip=19,stringsAsFactors=F)
 		names(g85)=c("g","reef_id","fit","phen","t","nmut","h","age","nad","mort")
		r=toupper(r)
		g85$reef_id=g85$reef_id+1
		g85$reef_id=as.factor(g85$reef_id)
		g85=subset(g85,g>=5300 & g<5900)
		g85=merge(g85,sst[,c(1,4)],by="reef_id",all.x=T)
		gg=g85 %>% mutate(cover=nad/km2)
		gg=gg %>% arrange(g)
		gg$dec=ceiling((gg$g-min(gg$g)+1)/10)
		cov.all=summarySE(gg, measurevar="cover",groupvars=c("reef_id","dec"))
		h.all=summarySE(gg, measurevar="h",groupvars=c("reef_id","dec"))
		age.all=summarySE(gg, measurevar="age",groupvars=c("reef_id","dec"))
		if (r=="A") { 
			covall=cov.all 
			ageall=age.all 
			hall=h.all
			} else { 
			  # adding empty age and cover values if the simulation differ in lengths
			  maxdec=max(cov.all$dec)
			    if (maxdec<max(covall$dec)){
  			    adds=subset(covall,dec>maxdec)
  			    adds$cover=0
  			    cov.all=rbind(cov.all,adds)
  			    cov.all=cov.all[order(cov.all$reef_id,cov.all$dec),]
  			    adds=subset(ageall,dec>maxdec)
  			    adds$age=0
  			    age.all=rbind(age.all,adds)
  			    age.all=age.all[order(age.all$reef_id,age.all$dec),]
  			    adds=subset(hall,dec>maxdec)
  			    adds$h=0
  			    h.all=rbind(h.all,adds)
  			    h.all=h.all[order(h.all$reef_id,h.all$dec),]
			    } 
  			  if (maxdec>max(covall$dec)){
  			    adds=subset(cov.all,dec>max(covall$dec))
			      adds$cover=0
			      covall=rbind(covall,adds)
			      covall=covall[order(covall$reef_id,covall$dec),]
			      adds=subset(age.all,dec>max(ageall$dec))
			      adds$age=0
			      ageall=rbind(ageall,adds)
			      ageall=ageall[order(ageall$reef_id,ageall$dec),]
			      adds=subset(h.all,dec>max(hall$dec))
			      adds$h=0
			      hall=rbind(hall,adds)
			      hall=hall[order(hall$reef_id,hall$dec),]
  			  }
  			covall[,4]=covall[,4]+cov.all[,4]
  			covall[,5]=covall[,5]+cov.all[,5]
  			ageall[,4]=ageall[,4]+age.all[,4]
  			ageall[,5]=ageall[,5]+age.all[,5]
  			hall[,4]=hall[,4]+h.all[,4]
  			hall[,5]=hall[,5]+h.all[,5]
			}
	}
	covall$dec=covall$dec-20	
	ageall$dec=ageall$dec-20	
	covall[,4]=covall[,4]/length(run)
	covall[,5]=covall[,5]/length(run)
	ageall[,4]=ageall[,4]/length(run)
	ageall[,5]=ageall[,5]/length(run)		
	save(hall,ageall, covall,file=paste(sets,".RData",sep=""))
}



