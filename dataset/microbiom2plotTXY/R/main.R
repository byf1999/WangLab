library(MicrobiomeAnalystR)
library(ggplot2)
library(readr)
#
check.rich<-function(input,step=100,enlarge=1){
  rich=as.data.frame(input)
  all_sampl=colnames(rich)
  reads_sum=apply(rich, 2, sum)
  size=max(reads_sum)*enlarge
  y=numeric()
  x=numeric()
  group=character()
  for (sampl in all_sampl) {
    all_reads=character()
    mappSet=rich[sampl]
    for (i in 1:nrow(mappSet)) {
      Gen=rownames(mappSet)[i]
      reads=rep(Gen,time=mappSet[i,1])
      all_reads=c(all_reads,reads)
    }
    Index=seq_along(all_reads)
    for (i in 0:ceiling(length(Index)/step)) {
      sele=sample(Index,size = i*step,replace = T)
      richness=length(unique(all_reads[sele]))
      x=c(x,step*i)
      y=c(y,richness)
      group=c(group,sampl)
    }
  }
  df=data.frame(X=x,Y=y,Group=group)
  library(ggplot2)
  p=ggplot(data = df,mapping = aes(x=X,y=Y,color=Group))+geom_smooth(se=F,formula = y~(log(x+1)))
  plot(p)
  return(df)
}
#
KEGG_abundance=function(mapper,targetLev,pathway_list,plot){
  mapper=read_file(mapper)
  pathway=read.delim(pathway_list,sep = '\t')
terms=unlist(strsplit(mapper,split = '\r\n'))
terms=terms[terms!='']
KO_index=grep('ko:',terms,value = F)
map_list=unlist(lapply(as.list(terms[-KO_index]), function(x){
  x=gsub('map','',x)
  return(x)
}))
Lev=c(rep('L3',each=nrow(pathway)),
      rep('L2',each=length(unique(pathway$L2))),
      rep('L1',each=length(unique(pathway$L1))))
names(Lev)=c(pathway$L3,unique(pathway$L2),unique(pathway$L1))
map_filtered=unlist(lapply(as.list(map_list), function(x){
  lis=unlist(strsplit(x,split = ' '))
  keyword=paste(lis[1:(length(lis)-1)],collapse = ' ')
  output=ifelse(Lev[keyword]%in%targetLev,paste(lis[1:(length(lis))],collapse = ' '),NA)
  return(output)
}))
map=na.omit(map_filtered)
abundance=unlist(lapply(as.list(map), function(x){
  lis=unlist(strsplit(x,split = ' '))
  output=as.numeric(gsub('[()]','',lis[length(lis)]))
  return(output)
}))
relative_abundance=abundance/mean(abundance)
tag=unlist(lapply(as.list(map), function(x){
  lis=unlist(strsplit(x,split = ' '))
  keyword=paste(lis[1:(length(lis)-1)],collapse = ' ')
  output=pathway$L1[which(pathway$L3==keyword)]
  return(output)
}))
df=data.frame(pathway=map,abundance=abundance,relative_abundance=relative_abundance,L1=tag)
df=df[order(df$relative_abundance,decreasing = T),]
if (plot){
  pdf=c()
  for (term in unique(df$L1)) {
    subdf=df[which(df$L1==term),]
    subdf$relative_abundance[subdf$relative_abundance>15]=NA
    subdf=na.omit(subdf)[1:5,]
    pdf=rbind(pdf,subdf)
  }
  p=ggplot(data=pdf,mapping = aes(relative_abundance,pathway))+
    geom_bar(stat = 'identity',aes(fill=L1))+theme_bw()+
    facet_grid(L1~.,scales = "free_y",space='free_y')+labs(fill='')
}else{p=NA}
output=list(df,p)
names(output)=c('data','plot')
return(output)
}
