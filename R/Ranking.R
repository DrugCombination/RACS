# This script is to performing manifold-ranking by calling LapRank function
Ranking <- function(DrugPairFeatureMatrix,Label,alpha){
  a <- DrugPairFeatureMatrix
  y <- Label
  # calculating drugpair similarity matrix using Eudian distance measure
  print('**calculating the drugpair similarity matrix...')
  w = as.matrix(dist(a,diag=T,upper=T))
  print('**Similarity matrix completed')
  print('Ranking...')
  r = LapRank_revised(w,y,alpha)
  print('Ranking completed')
  return(r)
}