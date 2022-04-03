# Created by Ruiyi Li on April 18, 2019
# This function return the value of objective function
#______________________
#Input:
# data: a numberic matrix with the rows are samples and colums are features, n_d *dim
# center: clustering centers with rows are cluster centers and colums are features, n_cl *dim
# cls: the cluster lable of samples, n_d *1
# weight: n_cl * dim
# delta: n_cl *1
# K: the number of clusters
# distType: COSINE, CITYBLOCK, euclidean

# Output: the value of objective function
#_______________________

objFun <- function(data,center,cls,weight,delta,K,distType)
{
  n_cl <- K  #the number of clusters
  n_d <- nrow(data) #the number of samples
  dim <- ncol(data) #the number of features
  
  temp <- weight*weight
  sum1 <- sum(delta[,1]*as.matrix(apply(temp,1,sum)))
  
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
    DW <-  weight1*abs(1.0/dim-data*center1)
  }else if(distType=="CITYBLOCK")
  {
    DW <-  weight1*abs(data-center1)
  }else
  {
    DW <-  weight1*(data-center1)^2
  }
  sum2 <- sum(apply(DW,1,sum))
  
  return(sum1+sum2)
}