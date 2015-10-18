# This is for testing of the impact of target information incomplete on model performance
# Used like query mode
# Function MI Defination===============================
MI.query <- function(queryDrugTargetVector,DrugTargetMatrix,GeneGOMatrix){
  query.DT <- queryDrugTargetVector
  DT.M <- DrugTargetMatrix
  GGO.M <- GeneGOMatrix
  DGO.Q <- query.DT %*% GGO.M
  DGO.Q[DGO.Q>0] = 1
  DGO.M <- DT.M %*% GGO.M
  DGO.M[DGO.M>0] = 1
  PX = rowSums(DGO.Q)/ncol(DGO.M)
  PY = rowSums(DGO.M)/ncol(DGO.M)
  PXY = (DGO.Q %*% t(DGO.M))/ncol(DGO.M)
  PXPY = PX %*% t(PY)
  M = PXY * log10(PXY/PXPY)
  M[is.na(M)] <- 0
  rownames(M) <- rownames(query.DT)
  colnames(M) <- rownames(DT.M)
  return(M)
}
# # MI demo test
# A <- matrix(c(1,0,0,0,1,1),nrow=2,byrow=T)
# rownames(A) <- c('D1','D2')
# colnames(A) <- c('T1','T2','T3')
# X <- t(as.matrix(A[1,]))
# rownames(X) <- rownames(A)[1]
# B <- matrix(c(0,1,0,0,0,1,1,1,0),nrow=3,byrow=T)
# MI.query(X,A,B)
#======================================================

# Dis==================================================
Dis.query <- function(queryDrugTargetVector,DrugTargetMatrix,InterTargetDisMatrix){
  query.DT <- queryDrugTargetVector
  DT.M <- DrugTargetMatrix
  ITD.M <- InterTargetDisMatrix
  N.drug <- nrow(DT.M)
  Dis.matrix <- matrix(0,nrow=nrow(query.DT),ncol=N.drug)
  rownames(Dis.matrix) <- rownames(query.DT)
  colnames(Dis.matrix) <- rownames(DT.M)
  for(i in 1:nrow(query.DT)){
    for(j in 1:N.drug){
      Dis.matrix[i,j] <- mean(InterTargetDisMatrix[query.DT[i,]>0,DT.M[j,]>0])
    }
  }
  return(Dis.matrix)
}
# # Dis demo test
# A <- matrix(c(0,1,0,1,0,1),nrow=2,byrow=T)
# rownames(A) <- c('D1','D2')
# colnames(A) <- c('T1','T2','T3')
# X <- t(as.matrix(A[1,]))
# rownames(X) <- rownames(A)[1]
# B <- matrix(c(0,2,3,2,0,1,3,1,0),nrow=3,byrow=T)
# Dis.query(X,A,B)
#======================================================
Efficiency <- function(g){
  sp <- shortest.paths(g)
  E <- mean(1/sp[upper.tri(sp)])
  return(E)
}
# DCI==================================================
DCI.query <- function(queryDrugTargetVector,DrugTargetMatrix,ppi.cancer,CN.nodes){
  query.DT <- queryDrugTargetVector
  DT.M <- DrugTargetMatrix
  Targets <- colnames(DT.M)
  Drugs <- rownames(DT.M)
  N.drug <- nrow(DT.M)
  E <- Efficiency(ppi.cancer)
  EX <- rep(0,nrow(query.DT))
  for(i in 1:nrow(query.DT)){
    target.x <- CN.nodes[CN.nodes[,2]%in%Targets[query.DT[i,]==1],1]
    gx <- ppi.cancer-target.x
    EX[i] <- Efficiency(gx)
  }
  #print(EX)
  EX.delta <- (E-EX)/E
  EY <- rep(0,N.drug)
  for(i in 1:N.drug){
    target.x <- CN.nodes[CN.nodes[,2]%in%Targets[DT.M[i,]==1],1]
    gx <- ppi.cancer-target.x
    EY[i] <- Efficiency(gx)
  }
  #print(EX)
  EY.delta <- (E-EY)/E
  EXY <- matrix(0,nrow=nrow(query.DT),ncol=N.drug)
  for(i in 1:nrow(query.DT)){
    target.x <- CN.nodes[CN.nodes[,2]%in%Targets[query.DT[i,]==1],1]
    print(i)
    for(j in 1:N.drug){
      target.y <- CN.nodes[CN.nodes[,2]%in%Targets[DT.M[j,]==1],1]
      target.xy <- union(target.x,target.y)
      gxy <- ppi.cancer - target.xy
      EXY[i,j] <- Efficiency(gxy)
    }
  }
  EXY.delta <- (E-EXY)/E
  DCI.matrix <- EXY.delta-as.matrix(EX.delta)%*%rep(1,length(EY.delta))-as.matrix(rep(1,length(EX.delta)))%*%EY.delta
  print(dim(DCI.matrix))
  print(length(Drugs))
  colnames(DCI.matrix) <- Drugs
  rownames(DCI.matrix) <- rownames(query.DT)
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
Eff.query <- function(queryDrugTargetVector,DrugTargetMatrix,ALL.nodes,CN.nodes,w=c(Degree.Global,Betweenness.Global,evcent.Global),lamda=0.1){
  query.DT <- queryDrugTargetVector
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
  Eff.w <- matrix(0,nrow=nrow(query.DT),ncol=N.drug)
  rownames(Eff.w) <- rownames(query.DT)
  colnames(Eff.w) <- rownames(DT.M)
  D <- sum(w)
  for(i in 1:nrow(query.DT)){
    for(j in 1:N.drug){
      TXY <- query.DT[i,] | DT.M[j,]
      CN <- TXY & t.CN
      BG <- TXY & t.BG
      NCN <- BG - CN
      A <- sum(t.w[CN])
      B <- sum(t.w[BG])
      C <- sum(t.w[NCN])
      Eff.w[i,j] <- lamda*(A/B)-(1-lamda)*(C/D)
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
MP.U.query <- function(queryDrugTargetVector,DrugTargetMatrix,TargetPathwayMatrix,PathwayUnRelatedMatrix){
  query.DT <- queryDrugTargetVector
  DT.M <- DrugTargetMatrix
  TP.M <- TargetPathwayMatrix
  PU.M <- PathwayUnRelatedMatrix[colnames(TP.M),colnames(TP.M)]
  DP.Q <- query.DT %*% TP.M
  DP.Q[DP.Q>0] <- 1
  DP.M <- DT.M %*% TP.M
  DP.M[DP.M>0] <- 1
  MPU <- DP.Q %*% (PU.M %*% t(DP.M)) / (rowSums(DP.Q) %*% t(rowSums(DP.M)))
  MPU[is.na(MPU)] <- 0
  return(MPU)
}
# print('***MP.U score calculating ...***')
# MP.U.matrix <- MP.U(DrugTargetMatrix,TargetPathwayMatrix,PathwayUnrelatedMatrix)
# print('***MP.U score complete!***')

##==================================================================================
##==================================================================================
ROCS.query <- function(queryDrugTargetVector,DrugTargetMatrix,GeneGOMatrix,InterTargetDisMatrix,PPI.Cancer,CN.nodes,
                 ALL.nodes,Degree.Global,Betweenness.Global,Evcent.Global,TargetPathwayMatrix,
                 PathwayUnrelatedMatrix){
  # calculate the score of each feature for all drug pairs
  print('**Step.1 calculate score of each feature for all drug pairs**')
  print('***MI score calculating ...***')
  MI.matrix <- MI.query(queryDrugTargetVector,DrugTargetMatrix,GeneGOMatrix)
  print('***MI score complete!***')
  print('***Dis score calculating ...***')
  Dis.matrix <- Dis.query(queryDrugTargetVector,DrugTargetMatrix,InterTargetDisMatrix)
  print('***Dis score complete!***')
  print('***DCI score calculating ...***')
  DCI.matrix <- DCI.query(queryDrugTargetVector,DrugTargetMatrix,PPI.Cancer,CN.nodes)
  print('***DCI score complete!***')
  print('***Eff score calculating ...***')
  print('***==Degree==***')
  Eff.D.matrix <- Eff.query(queryDrugTargetVector,DrugTargetMatrix,ALL.nodes,CN.nodes,w=Degree.Global)
  print('***==Degree Finish==***')
  print('***==Betweenness==***')
  Eff.B.matrix <- Eff.query(queryDrugTargetVector,DrugTargetMatrix,ALL.nodes,CN.nodes,w=Betweenness.Global)
  print('***==Betweenness Finish==***')
  print('***==Ev.cent.==***')
  Eff.E.matrix <- Eff.query(queryDrugTargetVector,DrugTargetMatrix,ALL.nodes,CN.nodes,w=Evcent.Global)
  print('***==Ev.cent. Finish==***')
  print('***Eff score complete!***')
  print('***MP.U score calculating ...***')
  MP.U.matrix <- MP.U.query(queryDrugTargetVector,DrugTargetMatrix,TargetPathwayMatrix,PathwayUnrelatedMatrix)
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
## query 2-target drugs target knock out
# query <- rbind(DrugTargetMatrix[rowSums(DrugTargetMatrix)==2,],DrugTargetMatrix[rowSums(DrugTargetMatrix)==2,])
# rownames(query) <- paste(rownames(query),c(rep(1,nrow(query)/2),rep(2,nrow(query)/2)),sep='.')
# for(i in 1:(nrow(query)/2))
#   query[i,query[i,]==1] <- c(0,1)
# query[(nrow(query)/2+1):nrow(query),]<- query[(nrow(query)/2+1):nrow(query),]-query[1:(nrow(query)/2),]
# if(!file.exists('data/queryTest'))
#   dir.create('data/queryTest')
# write.table(query,file='data/queryTest/query.2-targetdrug.txt',sep='\t',row.names=T,col.names=NA,quote=F)
######################################################
#query 3-target drugs target knock out
query.3 <- rbind(DrugTargetMatrix[rowSums(DrugTargetMatrix)==3,],DrugTargetMatrix[rowSums(DrugTargetMatrix)==3,],
                 DrugTargetMatrix[rowSums(DrugTargetMatrix)==3,],DrugTargetMatrix[rowSums(DrugTargetMatrix)==3,],
                 DrugTargetMatrix[rowSums(DrugTargetMatrix)==3,],DrugTargetMatrix[rowSums(DrugTargetMatrix)==3,])
rownames(query.3) <- paste(rownames(query.3),c(rep(1,nrow(query.3)/6),rep(2,nrow(query.3)/6),
                                               rep(3,nrow(query.3)/6),rep(4,nrow(query.3)/6),
                                               rep(5,nrow(query.3)/6),rep(6,nrow(query.3)/6)),sep='.')
N.query <- nrow(query.3)
for(i in 1:(N.query/6))
  query.3[i,query.3[i,]==1] <- c(0,1,1)
for(i in ((1:(N.query/6))+N.query/6*1))
  query.3[i,query.3[i,]==1] <- c(1,0,1)
for(i in ((1:(N.query/6))+N.query/6*2))
  query.3[i,query.3[i,]==1] <- c(1,1,0)
for(i in ((1:(N.query/6))+N.query/6*3))
  query.3[i,query.3[i,]==1] <- c(1,0,0)
for(i in ((1:(N.query/6))+N.query/6*4))
  query.3[i,query.3[i,]==1] <- c(0,1,0)
for(i in ((1:(N.query/6))+N.query/6*5))
  query.3[i,query.3[i,]==1] <- c(0,0,1)
if(!file.exists('data/queryTest'))
  dir.create('data/queryTest')
write.table(query.3,file='data/queryTest/query.3-targetdrug.txt',sep='\t',row.names=T,col.names=NA,quote=F)

#######
system.time({
  Score.matrix.3 <- ROCS.query(query.3[105:130,],DrugTargetMatrix,GeneGOMatrix,InterTargetDisMatrix,PPI.Cancer,
                       CN.nodes,ALL.nodes,Degree.Global,Betweenness.Global,
                       evcent.Global,TargetPathwayMatrix,PathwayUnrelatedMatrix)
})
SM.names <- names(Score.matrix.3)
for(i in 1:length(Score.matrix.3)){
  write.table(melt(t(Score.matrix.3[[i]])),file=paste("data/queryTest/5",SM.names[i],sep='/'),sep='\t',row.names=F,col.names=F,quote=F)
}

## evaluate the performance of RACS when target information is incomplete
# step 2. Ranking all drugpair
F1.test <- read.table('data/queryTest/MI.matrix',sep='\t',header=F)
F2.test <- read.table('data/queryTest/Dis.matrix',sep='\t',header=F)
F3.test <- read.table('data/queryTest/DCI.matrix',sep='\t',header=F)
F4.test <- read.table('data/queryTest/Eff.Degree.matrix',sep='\t',header=F)
F5.test <- read.table('data/queryTest/Eff.Betweenness.matrix',sep='\t',header=F)
F6.test <- read.table('data/queryTest/Eff.Evcent.matrix',sep='\t',header=F)
F7.test <- read.table('data/queryTest/MP.U.matrix',sep='\t',header=F)
DrugPairFeatureMatrix.test <- cbind(F1.test[,3],F2.test[,3],F3.test[,3],F4.test[,3],F5.test[,3],F6.test[,3],F7.test[,3])
DrugPairFeatureMatrix.test[is.infinite(DrugPairFeatureMatrix.test)] <- 0
write.table(DrugPairFeatureMatrix.test,file='data/queryTest/similarity.txt',sep='\t',row.names=F,col.names=F,quote=F)
colnames(DrugPairFeatureMatrix.test) <- c('MI','Dis','DCI','Eff.Degree','Eff.Betweenness',
                                     'Eff.Evcent','MP.U')
TEST.INDEX <- read.table('data/queryTest/TEST.INDEX.txt',sep='\t',header=F)
L1 <- read.table('data/26PositiveLabel.txt',sep='\t',header=F)
Label <- L1[,3]
alpha <- 0.9
N_DRUG <- nrow(DrugTargetMatrix)

for(i in 1:nrow(query.3)){
  A_I <- 1:(N_DRUG-1)+ (i-1)*(N_DRUG-1)
  DrugPairFeatureMatrix_i <- DrugPairFeatureMatrix
  DrugPairFeatureMatrix_i[TEST.INDEX[A_I,3],] <- DrugPairFeatureMatrix.test[A_I,]
  toSelect <- c(which(Label==1),TEST.INDEX[A_I,3])
  DPFM_i <- DrugPairFeatureMatrix_i[toSelect,]
  Rank.final2.test <- Ranking(DPFM_i,Label[toSelect],alpha)
  Rank.final.test <- cbind(L1[toSelect,],Rank.final2.test,DrugPairFeatureMatrix_i[toSelect,])
  write.table(Rank.final.test,file=paste('data/queryTest/',rownames(query.3)[i],'_00Result_Rank.txt',sep=''),sep='\t',row.names=F,col.names=T,quote=F)
  print(i)
}



