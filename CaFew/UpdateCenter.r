updateCenter <- function(data,cls,K,center)
{
  n_cl <- K
  dim <- ncol(data)
  n_d <- nrow(data)
  pre_center <- center
  
  center <- matrix(0,n_cl,dim)
  for (i in 1:n_cl)
  {
    loc <- which(cls==i)
    if(length(loc)>1)
    {
      tempD <- as.matrix(data[loc,])
      center[i,] <- apply(tempD,2,median)
    }else
    {
      center[i,] <- pre_center[i,]
      #locX <- floor(runif(1,1,n_d))
      #center[i,] <- data[locX,]
    }
  }
  
  
  return(center)
}