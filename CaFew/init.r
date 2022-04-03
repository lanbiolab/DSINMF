# Created by Ruiyi Li on April 18, 2019
# This function return the initialization value of center,delta,weight
#______________________
#Input:
# data: a numberic matrix with the rows are samples and colums are features, n_d *dim
# K : the number of clusters
# distType: COSINE, CITYBLOCK, euclidean
# k_delta

#Output
# center: clustering centers with rows are cluster centers and colums are features, n_cl *dim
# cls: the cluster lable of samples, n_d *1
# weight: n_cl * dim
# dropC: n_cl * dim
# delta: n_cl *1
#_______________________
find <- function(x)
{
  return(which(x==max(x))[1])
}

Init <- function(data,K,distType,k_delta)
{
   n_cl <- K
   dim <- ncol(data)
   n_d <- nrow(data)
   # initialize weight
   weight <- matrix(1/dim,n_cl,dim)
   #initialize center and cls
   #locX <- floor(runif(n_cl,1,n_d))
   #center  <- data[locX,]
   center <- InitCenter(data,n_cl)
   center <- as.matrix(center)
   #initialize cls
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
  
  # initialize delta
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
    delta1 <- weight1*abs(data-center1)
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
  
  
  iniValue <- list(center=center,cls=cls,weight=weight,delta=delta)
  return(iniValue)
}