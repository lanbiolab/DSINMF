find <- function(x)
{
  return(which(x==max(x))[1])
}

Patition <- function(data,center,weight,distType)
{
  n_cl <- nrow(center)
  dim <- ncol(data)
  n_d <- nrow(data)
  
  partition0 <- matrix(0,n_d,1)
  for (i in 1:n_cl)
  {
    center1 <- matrix(0,dim,n_d)
    weight1 <- matrix(0,dim,n_d)
    center1[,1:n_d] <- as.matrix(center[i,])
    center1 <- t(center1)
    weight1[,1:n_d] <- as.matrix(weight[i,])
    weight1 <- t(weight1)
    if(distType=="COSINE")
    {
      DW <- weight1*abs(1/dim-data*center1)
      temp <- as.matrix(apply(DW,1,sum))
      partition0 <- cbind(partition0,temp)
    }else if(distType=="CITYBLOCK")
    {
      DW <- weight1*abs(data-center1)
      temp <- as.matrix(apply(DW,1,sum))
      partition0 <- cbind(partition0,temp)
    }else
    {
      DW <- weight1*(abs(data-center1))^2
      temp <- as.matrix(apply(DW,1,sum))
      partition0 <- cbind(partition0,temp)
    }
  }
  partition0 <- partition0[,-1]
  partition0 <- as.matrix(partition0)
  cls <- as.matrix(apply(partition0,1,find))
  
  return(cls)
}