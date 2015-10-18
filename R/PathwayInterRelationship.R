# to build the pathway-pathway relation matrix
PathTargetMatrix <- as.matrix(read.table('data/Pathway-Target-Matrix.txt',sep='\t',header=T,row.names=1))
colnames(PathTargetMatrix) <- substr(colnames(PathTargetMatrix),1,8)
path.crosstalk <- t(PathTargetMatrix)%*%PathTargetMatrix
path.crosstalk[path.crosstalk>0] <- 1
hsa.targets <- rownames(PathTargetMatrix)
target.relationship <- read.csv('data/Gene-GeneInteraction.csv',header=F,colClasses='character')
targets.relation.matrix <- matrix(0,nrow=length(hsa.targets),ncol=length(hsa.targets))
rownames(targets.relation.matrix) <- hsa.targets
colnames(targets.relation.matrix) <- hsa.targets
for(i in 1:nrow(target.relationship)){
  if(target.relationship[i,1]%in%hsa.targets & target.relationship[i,2]%in%hsa.targets){
    targets.relation.matrix[target.relationship[i,1],target.relationship[i,2]] <- 1
    targets.relation.matrix[target.relationship[i,2],target.relationship[i,1]] <- 1
  }
}
path.interacting <- t(PathTargetMatrix) %*% t(t(PathTargetMatrix) %*% targets.relation.matrix)
diag(path.interacting) <- 0
path.interacting[path.interacting>0] <- 1
path.interacting[path.crosstalk>0] <- 0

path.Unrelated <- path.crosstalk + path.interacting
path.Unrelated[path.Unrelated>0] <- -1
path.Unrelated[path.Unrelated==0] <- 1
path.Unrelated[path.Unrelated<0] <- 0

write.table(path.Unrelated,file='data/path.Unrelated.matrix.txt',sep='\t',row.names=T,col.names=NA,quote=F)
write.table(path.crosstalk,file='data/path.crosstalk.matrix.txt',sep='\t',row.names=T,col.names=NA,quote=F)
write.table(path.interacting,file='data/path.interacting.matrix.txt',sep='\t',row.names=T,col.names=NA,quote=F)

