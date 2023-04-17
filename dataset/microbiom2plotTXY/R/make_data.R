
L1map=unlist(strsplit(read_file('L1pathway.txt'),split = '\r\n'))

L2map=unlist(strsplit(read_file('L2+L1pathway.txt'),split = '\r\n'))
L2map=L2map[which(!L2map%in%L1map)]

allmap=unlist(strsplit(read_file('L3+L2+L1 pathway.txt'),split = '\r\n'))
allmap=unique(allmap)
allmap=allmap[allmap!='']
L3map=allmap[which(!allmap%in%c(L1map,L2map))]
L1index=which(allmap%in%L1map)
L2index=which(allmap%in%L2map)
L1=c()
L2=c()
for (map in L3map) {
  i=which(allmap==map)
  l2_i=sort(L2index[L2index<i],decreasing = T)[1]
  l1_i=sort(L1index[L1index<i],decreasing = T)[1]
  L1=c(L1,allmap[l1_i])
  L2=c(L2,allmap[l2_i])
}



df=as.data.frame(cbind(L3=L3map,L2,L1))
write_delim(df,'pathway_list.txt',delim = '\t')

#




