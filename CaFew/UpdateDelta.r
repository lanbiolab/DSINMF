# Created by Ruiyi Li on April 18, 2019
# This function return the value of updated weight
#______________________
#Input:
# data: a numberic matrix with the rows are samples and colums are features, n_d *dim
# pre_center: clustering centers with rows are cluster centers and colums are features, n_cl *dim
# pre_cls: the cluster lable of samples, n_d *1
# pre_weight: n_cl * dim
# K: the number of clusters
# distType: COSINE, CITYBLOCK, euclidean
# k_delta

# Output: the updated delta
#_______________________


updateDelta <- function(data,center,cls,weight,K,k_delta)
{
  n_cl <- K
  dim <- ncol(data)
  n_d <- nrow(data)
  
  
  ss_v <- as.matrix(apply(weight^2,1,sum))
  
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
    delta1 <- weight1*abs(1.0/dim-data*center1)
  }else if(distType=="CITYBLOCK")
  {
    delta1 <- weight1* abs(data-center1)
  }else
  {
    delta1 <- weight1*(data-center1)^2
  }
  delta1 <- as.matrix(apply(delta1,1,sum))
  
  delta <- matrix(0,n_cl,1)
  for (i in 1:n_cl)
  {
    loc <- which(cls[,1]==i)
    if(length(loc)>0)
    {
      delta[i,1] <- sum(delta1[loc,],na.rm=T)
    }
  }
  for(i in 1:n_cl)
  {
    if(ss_v[i,1]!=0)
    {
      delta[i,1] <- k_delta*delta[i,1]/ss_v[i,1] 
    }
  }
  
  return(delta)
}