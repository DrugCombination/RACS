# To build target-target distance matrix
library(igraph)
targets <- as.vector(read.table('data/all_target.txt',sep='\t')[,1])
PPI <- read.graph('data/6_ppi.net',format='pajek')
nodes <- read.table('data/6_ppi.node',sep=' ',header=F)
targets.index <- rep(0,length(targets))
for(i in 1:length(targets)){
  if(targets[i] %in% nodes[,2])
    targets.index[i] <- nodes[nodes[,2]==targets[i],1]
}
Shortest.path <- shortest.paths(PPI,v=targets.index[targets.index!=0],to=targets.index[targets.index!=0])
rownames(Shortest.path) <- targets[targets.index!=0]
colnames(Shortest.path) <- targets[targets.index!=0]
SP <- matrix(0,nrow=length(targets),ncol=length(targets))
rownames(SP) <- targets
colnames(SP) <- targets
for(i in rownames(SP)){
  for(j in rownames(SP)){
    if(i %in% rownames(Shortest.path) && j %in% rownames(Shortest.path))
      SP[which(rownames(SP)==i),which(rownames(SP)==j)] <- Shortest.path[which(rownames(Shortest.path)==i),which(rownames(Shortest.path)==j)]
    else
      SP[which(rownames(SP)==i),which(rownames(SP)==j)] <- 0
  }
}
write.table(SP,file='data/inter-targetShortestPathLength.matrix2.txt',sep='\t',row.names=T,col.names=T,quote=F)
