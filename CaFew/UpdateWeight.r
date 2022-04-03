# Created by Ruiyi Li on April 18, 2019
# This function return the value of updated weight
#______________________
#Input:
# data: a numberic matrix with the rows are samples and colums are features, n_d *dim
# center: clustering centers with rows are cluster centers and colums are features, n_cl *dim
# cls: the cluster lable of samples, n_d *1
# weight: n_cl * dim
# delta: n_cl *1
# K: the number of clusters
# distType: COSINE, CITYBLOCK, euclidean
#dropC: n_cl * dim, drop-out rate of genes in different clusters

# Output: the updated weight
#_______________________



updateWeight <- function(weight,data,center,cls,delta,K,distType)
{
  
  n_cl <- K
  dim <- ncol(data)
  n_d <- nrow(data)
  
  
  center1 <- matrix(0,dim,n_d)
  weight1 <- matrix(0,dim,n_d)
  for (i in 1:n_cl)
  {
    loc <- which(cls[,1]==i)
    if(length(loc)>0)
    {
      center1[,loc] <- as.matrix(center[i,])
      weight1[,loc] <- as.matrix(weight[i,])
    }
  }
  center1 <- t(center1)
  weight1 <- t(weight1)
  
  if(distType=="COSINE")
  {
    DW <- abs(1.0/dim-data*center1)
  }else if(distType=="CITYBLOCK")
  {
    DW <- abs(data-center1)
  }else
  {
    DW <- (data-center1)^2
  }
  
  DWT <- DW*weight1
  consD <- as.matrix(apply(DWT,1,sum))
  consD <- consD/dim
  DW <- consD[,1]-DW
  
  weight1 <- matrix(0,n_cl,dim)
  for (i in 1:n_cl)
  {
    loc <- which(cls[,1]==i)
    if(length(loc)>0)
    {
      weight1[i,] <- apply(as.matrix(DW[loc,]),2,sum)
    }
  }
  delta[which(delta[,1]==0),] <- 0.001
  deltaT <- matrix(0,n_cl,dim)
  deltaT[,1:dim] <- (1/(2*delta))[,1]
  weight1 <- 1/dim + deltaT*weight1
  
  
  #modify negative weights
  for (i in 1:n_cl)
  {
    minX <- min(weight1[i,])
    if(minX<0)
    {
      weight1[i,which(weight1[i,]<0)] <- weight1[i,which(weight1[i,]<0)]-minX
    }
  }
  sum_v <- as.matrix(apply(weight1,1,sum))
  for(i in 1:n_cl)
  {
    if(sum_v[i,1]!=0)
    {
      weight1[i,] <- weight1[i,]/sum_v[i,1]
    }
  }
  
  weight <- weight1
  
  temp <- weight*weight
  sum1 <- sum(as.matrix(apply(temp,1,sum)))
  
  return(list(weight=weight,he=sum1))
}