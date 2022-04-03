source("InitCenter.r")
source("init.r")
source("updateWeight.r")
source("Patition.r")
source("updateCenter.r")
source("updateDelta.r")
source("objFun.r")

find0 <- function(x)
{
  return(length(which(x==0))/length(x))
}
findL <- function(x)
{
  loc <- which(x>=0.99)
  temp=0
  if (length(loc)>0)
  {
    temp=1
  }
  return(temp)
}

data_name = 'mouse2'

data <- read.table(file = paste0("E:/single_data/",data_name,"_ready.txt"))
trueclass <- read.table(file = paste0("E:/single_data/",data_name,"_ready_label.txt"))
k = length(unique(trueclass))

c0cnt <- apply(data,2,find0)
loc <- which(c0cnt>0.98)
if(length(loc)>0)
{
  data <- data[,-loc]
}

re <- cor(data)
re[row(re)>col(re)] <- 0 
diag(re) <- 0
gn <- nrow(re)
re1 <- apply(re,1,findL)
loc <- which(re1==1)
if(length(loc)>0)
{
  data <- data[,-loc]
}

distType="CITYBLOCK"
k_delta=0.001
inRe <- Init(data,K,distType,k_delta)
center <- inRe$center
cls <- inRe$cls
weight <- inRe$weight
delta <- inRe$delta
iterR1 <- matrix(0,50,2)
i <- 1
while (i < 51)
{
  re <- updateWeight(weight,data,center,cls,delta,K,distType)
  iterR1[i,2] <- re$he
  weight <- re$weight
  cls <- Patition(data,center,weight,distType)
  center <- updateCenter(data,cls,K,center)
  delta <- updateDelta(data,center,cls,weight,K,k_delta)
  iterR1[i,1] <- objFun(data,center,cls,weight,delta,K,distType)
  i <- i+1
}

weight <- as.matrix(weight)
dim <- ncol(weight)
max_w <- as.matrix(apply(weight,2,max))
hr <- hist(max_w)
th1 <- hr$breaks[2]
RetainIndex1 <- which(max_w[,1]>=th1)
data1 <- data[,RetainIndex1] # The result after the first selecting step


weight1 <- weight[,RetainIndex1]
while(ncol(data1)>10000)
{
  max_w <- max_w[RetainIndex1]
  hr <- hist(max_w)
  th1 <- hr$breaks[2]
  RetainIndex1 <- which(max_w>=th1)
  data1 <- data1[,RetainIndex1]
  weight1 <- weight1[,RetainIndex1]
}


y <- apply(weight1,2,sd)
x <- apply(weight1,2,mean)
CV <- y/x 
fR <- lm(log(CV^2)~log10(x))
difV <- fR$residuals
p_value <- pnorm(-abs(scale(difV, center=TRUE,scale=TRUE)))
RetainIndex2 <- which(p_value<=0.05)
data2 <- data1[,RetainIndex2] # The result after the second selecting step
write.csv(data2, paste0("E:/single_data/featureSelection/",data_name,"_data_p_value005.csv"))
write.csv(p_value, paste0("E:/single_data/featureSelection/",data_name,"_p_value.csv"))
p <- order(p_value)
write.csv(p, paste0("E:/single_data/featureSelection/",data_name,"_p_value_order.csv"))


write.csv(data2, paste0("E:/single_data/featureSelection/mouse2_data_p_value005.csv"))
write.csv(p_value, paste0("E:/single_data/featureSelection/mouse2_p_value.csv"))
p <- order(p_value)
write.csv(p, paste0("E:/single_data/featureSelection/mouse2_p_value_order.csv"))

