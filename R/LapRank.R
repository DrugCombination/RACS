
LapRank_revised <- function(W, y, alpha){
  ###
  # Laplacian-based Ranking
  # Zhen, Sheng: shengzhen@live.cn
  # INPUT: 
  #   (1) W: n*n, similarity matrix
  #   (2) y: n*1, label vector: 1:positive, 0: unknown
  #   (3) alpha: tradeoff parameter --> tradeoff_alpha = alpha/n
  # OUTPUT: 
  #   (1) f: n*1, score of each instance
  ###
  n = length(y)
  tradeoff_alpha = alpha
  # --- Normalized Laplacian Matrix
  print('step.1')
  D <- diag(rowSums(W))
  print('step.2')
  DA <- D^(-0.5)
  print('step.3')
  DA[is.infinite(DA)] <- 0
  print('step.4')
  L <- DA %*% W %*% DA
  print('step.5')
  FA = diag(n) - tradeoff_alpha * L
  print('step.6')
  FB = (1-tradeoff_alpha) * y
  print('step.7')
  f = solve(FA,FB)
  return(f)
}
