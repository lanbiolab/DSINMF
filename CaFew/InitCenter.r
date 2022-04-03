# library("ClusterX")
InitCenter <- function(data,K)
{
  n_cl <- K
  dim <- ncol(data)
  n_d <- nrow(data)
  
  # ans <- ClusterX(data,outDim = 50)
  ans <- kmeans(data,n_cl)
  #rho <- (ans$rho)/max(ans$rho)
  #delta <- (ans$delta)/max(ans$delta)
  #re <- as.matrix(rho*delta)
  #loc <- order(re,decreasing = T)[1:K]
  #center <- data[loc,]
  return(ans$centers)
}