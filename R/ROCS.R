# Function MI Defination===============================
MI <- function(DrugTargetMatrix,GeneGOMatrix){
  DT.M <- DrugTargetMatrix
  GGO.M <- GeneGOMatrix
  DGO.M <- DT.M %*% GGO.M
  DGO.M[DGO.M>0] = 1
  PX = rowSums(DGO.M)/ncol(DGO.M)
  PXY = (DGO.M %*% t(DGO.M))/ncol(DGO.M)
  PXPY = PX %*% t(PX)
  M = PXY * log10(PXY/PXPY)
  rownames(M) <- rownames(DT.M)
  colnames(M) <- rownames(DT.M)
  M[is.na(M)] <- 0
  return(M)
}
# MI demo test
# A <- matrix(c(1,0,0,0,1,1),nrow=2,byrow=T)
# B <- matrix(c(0,1,0,1,0,0,1,0,1,0,0,1),nrow=3,byrow=T)
# MI(A,B)
#======================================================

# Dis==================================================
Dis <- function(DrugTargetMatrix,InterTargetDisMatrix){
  DT.M <- DrugTargetMatrix
  ITD.M <- InterTargetDisMatrix
  N.drug <- nrow(DT.M)
  Dis.matrix <- matrix(0,nrow=N.drug,ncol=N.drug)
  rownames(Dis.matrix) <- rownames(DT.M)
  colnames(Dis.matrix) <- rownames(DT.M)
  for(i in 1:N.drug){
    for(j in i:N.drug){
      Dis.matrix[i,j] <- mean(InterTargetDisMatrix[DT.M[i,]>0,DT.M[j,]>0])
      Dis.matrix[j,i] <- Dis.matrix[i,j]
    }
  }
  return(Dis.matrix)
}
# Dis demo test
# A <- matrix(c(0,1,0,1,0,1),nrow=2,byrow=T)
# B <- matrix(c(0,2,3,2,0,1,3,1,0),nrow=3,byrow=T)
# Dis(A,B)
#======================================================
Efficiency <- function(g){
  sp <- shortest.paths(g)
  E <- mean(1/sp[upper.tri(sp)])
  return(E)
}
# DCI==================================================
DCI <- function(DrugTargetMatrix,ppi.cancer,CN.nodes){
  DT.M <- DrugTargetMatrix
  Targets <- colnames(DT.M)
  Drugs <- rownames(DT.M)
  N.drug <- nrow(DT.M)
  E <- Efficiency(ppi.cancer)
  EX <- rep(0,N.drug)
  for(i in 1:N.drug){
    target.x <- CN.nodes[CN.nodes[,2]%in%Targets[DT.M[i,]==1],1]
    gx <- ppi.cancer-target.x
    EX[i] <- Efficiency(gx)
  }
  #print(EX)
  EX.delta <- (E-EX)/E
  EXY <- matrix(0,nrow=N.drug,ncol=N.drug)
  for(i in 1:N.drug){
    target.x <- CN.nodes[CN.nodes[,2]%in%Targets[DT.M[i,]==1],1]
    print(i)
    for(j in i:N.drug){
      target.y <- CN.nodes[CN.nodes[,2]%in%Targets[DT.M[j,]==1],1]
      target.xy <- union(target.x,target.y)
      gxy <- ppi.cancer - target.xy
      EXY[i,j] <- Efficiency(gxy)
      EXY[j,i] <- EXY[i,j]
    }
  }
  EXY.delta <- (E-EXY)/E
  DCI.matrix <- EXY.delta-as.matrix(EX.delta)%*%rep(1,length(EX.delta))-as.matrix(rep(1,length(EX.delta)))%*%EX.delta
  colnames(DCI.matrix) <- Drugs
  rownames(DCI.matrix) <- Drugs
  return(DCI.matrix)
}
#system.time(XX <- DCI(DrugTargetMatrix,PPI.Cancer,CN.nodes))

#======================================================

# Efficacy (Degree, Betweenness, eigenvector centrality)
# Degree.Global <- degree(PPI)
# Betweenness.Global <- betweenness(PPI)
# evcent.Global <- evcent(PPI)$vector
# ALL.nodes # all proteins in backgroud PPI
# CN.nodes # all proteins in Cancer Network
Eff <- function(DrugTargetMatrix,ALL.nodes,CN.nodes,w=c(Degree.Global,Betweenness.Global,evcent.Global),lamda=0.1){
  DT.M <- DrugTargetMatrix
  Targets <- colnames(DT.M)
  N.drug <- nrow(DT.M)
  t.BG <- Targets %in% ALL.nodes[,2]
  t.CN <- Targets %in% CN.nodes[,2]
  t.w <- rep(0,length(Targets))
  for(i in 1:length(Targets)){
    if(Targets[i] %in% ALL.nodes[,2])
      t.w[i] <- w[ALL.nodes[ALL.nodes[,2]==Targets[i],1]]
  }
  Eff.w <- matrix(0,nrow=N.drug,ncol=N.drug)
  rownames(Eff.w) <- rownames(DT.M)
  colnames(Eff.w) <- rownames(DT.M)
  D <- sum(w)
  for(i in 1:N.drug){
    for(j in i:N.drug){
      TXY <- DT.M[i,] | DT.M[j,]
      CN <- TXY & t.CN
      BG <- TXY & t.BG
      NCN <- BG - CN
      A <- sum(t.w[CN])
      B <- sum(t.w[BG])
      C <- sum(t.w[NCN])
      Eff.w[i,j] <- lamda*(A/B)-(1-lamda)*(C/D)
      Eff.w[j,i] <- Eff.w[i,j]
    }
  }
  return(Eff.w)
}
# print('***Eff score calculating ...***')
# print('***==Degree==***')
# Eff.D.matrix <- Eff(DrugTargetMatrix,ALL.nodes,CN.nodes,w=Degree.Global)
# print('***==Degree Finish==***')
# print('***==Betweenness==***')
# Eff.B.matrix <- Eff(DrugTargetMatrix,ALL.nodes,CN.nodes,w=Betweenness.Global)
# print('***==Betweenness Finish==***')
# print('***==Ev.cent.==***')
# Eff.E.matrix <- Eff(DrugTargetMatrix,ALL.nodes,CN.nodes,w=evcent.Global)
# print('***==Ev.cent. Finish==***')
# print('***Eff score complete!***')
#===============================================

#MP.U
# pathway.list <- read.table('data/Pathway.txt',sep='\t',header=T,quote='\"')[,2]
# #pathway.UnRelated <- reat.table('data/path.Unrelated.matrix.txt',sep='\t',header=T,row.names=1)
# pathway.UnRelated <- path.Unrelated[pathway.list,pathway.list]
# protein_pathway.matrix <- read.table('data/Pathway-Target-Matrix.uniprot.txt',sep='\t',header=T,row.names=1)
# target_pathway.matrix <- matrix(0,nrow=length(targets),ncol=nrow(pathway.UnRelated))
# colnames(target_pathway.matrix) <- colnames(pathway.UnRelated)
# rownames(target_pathway.matrix) <- t(targets)
# for(i in 1:nrow(target_pathway.matrix)){
#   if(rownames(target_pathway.matrix)[i] %in% rownames(protein_pathway.matrix)){
#     for(j in 1:ncol(target_pathway.matrix)){
#     target_pathway.matrix[i,j] <- protein_pathway.matrix[rownames(target_pathway.matrix)[i],colnames(pathway.UnRelated)[j]]
#     }
#   }
# }
# write.table(target_pathway.matrix,file='data/target_pathway.matrix',sep='\t',row.names=T,col.names=NA,quote=F)
MP.U <- function(DrugTargetMatrix,TargetPathwayMatrix,PathwayUnRelatedMatrix){
  DT.M <- DrugTargetMatrix
  TP.M <- TargetPathwayMatrix
  PU.M <- PathwayUnRelatedMatrix[colnames(TP.M),colnames(TP.M)]
  DP.M <- DT.M %*% TP.M
  DP.M[DP.M>0] <- 1
  MPU <- DP.M %*% (PU.M %*% t(DP.M)) / (rowSums(DP.M) %*% t(rowSums(DP.M)))
  MPU[is.na(MPU)] <- 0
  return(MPU)
}
# print('***MP.U score calculating ...***')
# MP.U.matrix <- MP.U(DrugTargetMatrix,TargetPathwayMatrix,PathwayUnrelatedMatrix)
# print('***MP.U score complete!***')

##==================================================================================
##==================================================================================
ROCS <- function(DrugTargetMatrix,GeneGOMatrix,InterTargetDisMatrix,PPI.Cancer,CN.nodes,
                 ALL.nodes,Degree.Global,Betweenness.Global,Evcent.Global,TargetPathwayMatrix,
                 PathwayUnrelatedMatrix){
  # calculate the score of each feature for all drug pairs
  print('**Step.1 calculate score of each feature for all drug pairs**')
  print('***MI score calculating ...***')
  MI.matrix <- MI(DrugTargetMatrix,GeneGOMatrix)
  print('***MI score complete!***')
  print('***Dis score calculating ...***')
  Dis.matrix <- Dis(DrugTargetMatrix,InterTargetDisMatrix)
  print('***Dis score complete!***')
  print('***DCI score calculating ...***')
  DCI.matrix <- DCI(DrugTargetMatrix,PPI.Cancer,CN.nodes)
  print('***DCI score complete!***')
  print('***Eff score calculating ...***')
  print('***==Degree==***')
  Eff.D.matrix <- Eff(DrugTargetMatrix,ALL.nodes,CN.nodes,w=Degree.Global)
  print('***==Degree Finish==***')
  print('***==Betweenness==***')
  Eff.B.matrix <- Eff(DrugTargetMatrix,ALL.nodes,CN.nodes,w=Betweenness.Global)
  print('***==Betweenness Finish==***')
  print('***==Ev.cent.==***')
  Eff.E.matrix <- Eff(DrugTargetMatrix,ALL.nodes,CN.nodes,w=Evcent.Global)
  print('***==Ev.cent. Finish==***')
  print('***Eff score complete!***')
  print('***MP.U score calculating ...***')
  MP.U.matrix <- MP.U(DrugTargetMatrix,TargetPathwayMatrix,PathwayUnrelatedMatrix)
  print('***MP.U score complete!***')
  print('**Step.1 complete!**')
  SCORE <- list(MI.matrix=MI.matrix,Dis.matrix=Dis.matrix,DCI.matrix=DCI.matrix,
                Eff.Degree.matrix=Eff.D.matrix,
                Eff.Betweenness.matrix=Eff.B.matrix,
                Eff.Evcent.matrix=Eff.E.matrix,MP.U.matrix=MP.U.matrix)
  return(SCORE)
  print('--*--*--*--*--*--*--*--')
}
#================================================================================#
#================================================================================#
#================================================================================#
#================================================================================#
# test
library(igraph)
library(reshape2)
DrugTargetMatrix <- as.matrix(read.table('data/drug-target.matrix.txt',sep='\t',header=T,row.names=1))
GeneGOMatrix <- as.matrix(read.table('data/GENE-GO.matrix.txt',sep='\t',header=T,row.names=1))
InterTargetDisMatrix <- as.matrix(read.table('data/inter-targetShortestPathLength.matrix2.txt',sep='\t',header=T,row.names=1))
PPI.Cancer <- read.graph('data/Cancer_Network.net',format='pajek')
CN.nodes <- read.table('data/CN.node',sep=' ',header=F)
ALL.nodes <- read.table('data/6_ppi.node',sep=' ',header=F)
PPI <- read.graph('data/6_ppi.net',format='pajek')
Degree.Global <- degree(PPI)
Betweenness.Global <- betweenness(PPI)
evcent.Global <- evcent(PPI)$vector
TargetPathwayMatrix <- as.matrix(read.table('data/target_pathway.matrix',sep='\t',header=T,row.names=1))
PathwayUnrelatedMatrix <- as.matrix(read.table('data/path.Unrelated.matrix.txt',sep='\t',header=T,row.names=1))
system.time({
Score.matrix <- ROCS(DrugTargetMatrix,GeneGOMatrix,InterTargetDisMatrix,PPI.Cancer,
                     CN.nodes,ALL.nodes,Degree.Global,Betweenness.Global,
                     evcent.Global,TargetPathwayMatrix,PathwayUnrelatedMatrix)
})
SM.names <- names(Score.matrix)
for(i in 1:length(Score.matrix)){
  write.table(melt(Score.matrix[[i]])[upper.tri(Score.matrix[[i]]),],file=paste("data",SM.names[i],sep='/'),sep='\t',row.names=F,col.names=F,quote=F)
}
# step 2. Ranking all drugpair
F1 <- read.table('data/MI.matrix',sep='\t',header=F)
F2 <- read.table('data/Dis.matrix',sep='\t',header=F)
F3 <- read.table('data/DCI.matrix',sep='\t',header=F)
F4 <- read.table('data/Eff.Degree.matrix',sep='\t',header=F)
F5 <- read.table('data/Eff.Betweenness.matrix',sep='\t',header=F)
F6 <- read.table('data/Eff.Evcent.matrix',sep='\t',header=F)
F7 <- read.table('data/MP.U.matrix',sep='\t',header=F)
DrugPairFeatureMatrix <- cbind(F1[,3],F2[,3],F3[,3],F4[,3],F5[,3],F6[,3],F7[,3])
DrugPairFeatureMatrix[is.infinite(DrugPairFeatureMatrix)] <- 0
write.table(DrugPairFeatureMatrix,file='data/similarity.txt',sep='\t',row.names=F,col.names=F,quote=F)
colnames(DrugPairFeatureMatrix) <- c('MI','Dis','DCI','Eff.Degree','Eff.Betweenness',
                                     'Eff.Evcent','MP.U')
L1 <- read.table('data/26PositiveLabel.txt',sep='\t',header=F)
Label <- L1[,3]
alpha <- 0.9
Rank.final2 <- Ranking(DrugPairFeatureMatrix,Label,alpha)
Rank.final <- cbind(L1,Rank.final2,DrugPairFeatureMatrix)
write.table(Rank.final,file='data/00Result_Rank.txt',sep='\t',row.names=F,col.names=T,quote=F)







