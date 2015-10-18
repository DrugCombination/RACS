DrugTargetMatrix <- as.matrix(read.table('data/drug-target.matrix.txt',sep='\t',header=T,row.names=1))
pathway.list <- read.table('data/Pathway.txt',sep='\t',header=T,quote='\"')[,2]
pathway.UnRelated <- read.table('data/path.Unrelated.matrix.txt',sep='\t',header=T,row.names=1)
pathway.UnRelated <- pathway.UnRelated[pathway.list,pathway.list]
protein_pathway.matrix <- read.table('data/Pathway-Target-Matrix.uniprot.txt',sep='\t',header=T,row.names=1)
targets <- colnames(DrugTargetMatrix)
target_pathway.matrix <- matrix(0,nrow=length(targets),ncol=nrow(pathway.UnRelated))
colnames(target_pathway.matrix) <- colnames(pathway.UnRelated)
rownames(target_pathway.matrix) <- t(targets)
for(i in 1:nrow(target_pathway.matrix)){
  if(rownames(target_pathway.matrix)[i] %in% rownames(protein_pathway.matrix)){
    for(j in 1:ncol(target_pathway.matrix)){
      target_pathway.matrix[i,j] <- protein_pathway.matrix[rownames(target_pathway.matrix)[i],colnames(pathway.UnRelated)[j]]
    }
  }
}
write.table(target_pathway.matrix,file='data/target_pathway.matrix',sep='\t',row.names=T,col.names=NA,quote=F)
