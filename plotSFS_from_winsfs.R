sfs2matrix=function(sfs,n1,n2,invariants=FALSE) {
  dd=matrix(ncol=2*n1+1,nrow=2*n2+1,sfs)
  if (invariants==FALSE) { dd[1,1]=dd[2*n2+1,2*n1+1]=0 }
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

