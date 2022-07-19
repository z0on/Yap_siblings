pops=seq(1,3) # population names
e.local=c(-2,0,2) # local environmental settting
sin.p=5 # period of sinusoidal climate fluctuations, in reproduction cycles
# ---- rand /  sin settings:
# e.sd= 0.2 # SD of random environmentsl fluctuations (different in each pop)
sin.amp=0.5 # amplitude of sinusoidal climate fluctuations, common to all pops
# sin.rand=0.2 # amplitude of random climate fluctuations, common to all pops
e.sd= 0.0 # SD of random environmentsl fluctuations (different in each pop)
sin.rand=0.0 # amplitude of random climate fluctuations, common to all pops
burnin=5500 # length of burnin period
gmax=6000 # maximum number of generations

e.increment=2.38/90 # increment in environmental setting per reproduction cycle (year) past burn-in period

# for no-fluctuation and sinusoidal model, there is only one environmental profile ("a") for all SLiM replicates; for random model there are four different profiles.
for (index in c("a")) {
  message(index)
  dd=lapply(seq_len(gmax),function(i) {
    newgen=e.local+sin(i*2*pi/sin.p)*0.5*(sin.amp+rnorm(1,0,sin.rand))
    newgen=newgen+rnorm(length(newgen),0,e.sd)
    if (i>burnin){
      newgen=newgen+(i-burnin)*e.increment
    }
    return(round(newgen,3))
  }
  )
  envs=data.frame(do.call(cbind,dd))
  write.table(envs,row.names=F,quote=F,sep="\t",file=paste("toy_",index,"_environment.txt",sep=""))
}
ee=data.frame(t(envs[1:3,]))
plot(ee[,1],type="l",xlim=c(5400,6000))
lines(ee[,2],type="l",col="red")

