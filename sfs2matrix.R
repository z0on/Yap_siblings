
sfs2matrix=function(sfs,n1,n2) {
  dd=matrix(ncol=2*n1+1,nrow=2*n2+1,sfs)
#  dd[1,1]=dd[2*n2+1,2*n1+1]=0
  return(apply(dd,2,rev))
}


# sfs=scan("~/Dropbox/ok_36_26.sfs")
# n1=36
# n2=26
# sfs2d=sfs2matrix(sfs,n1,n2)
# pheatmap::pheatmap(log(sfs2d+0.5,10),cluster_rows = F, cluster_cols = F,border_color = F)
