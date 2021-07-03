# Tsallis Entropy based feature selection
Tsallisfeature <- function(data,cell,p,q,nf){
 library(foreach)
library(doParallel)
 
 #data=a matrix format, Cells sholud be in rwo, Genes should be in coloumn in data, header=False
 # cell = cell types of data, header=False
data=as.matrix(data)
classc=cell
n <- nrow(data)
col<-ncol(data)
count=ncol(data)
#nf=Number of feature to be selected.
#p=Number of cores. 
#q= Tuned parameter of Tsallis.
 
 # Tsallis Entropy
cl <- makeCluster(p)
registerDoParallel(cl)
fea<- matrix(0, nrow=1,ncol =nf)
parl<-foreach(j=1:count, .combine=c) %dopar%
  {
    u<-as.matrix(cbind(classc,data[,j]))
    mytable=table(u)
    s1=as.data.frame(prop.table(mytable))
    nem=sum(s1$Freq^q)
    mytable=table(u[,1])
    s1=as.data.frame(prop.table(mytable))
    denom=sum(s1$Freq^q)
    trans=(1/(q-1))*(1-(nem/denom))
  }

mimc1<-as.matrix(parl)
fea[1,1]=which.max(mimc1)
idx=which.max(mimc1)
for(i in 1:(nf-1))
{  
  parl<-foreach(j=1:count, .combine=c) %dopar%
    {
      u2<-as.matrix(cbind(classc,data[,idx],data[,j]))
      mytable=table(u2)
      s1=as.data.frame(prop.table(mytable))
      nem=sum(s1$Freq^q)
      mytable=table(u2[,2],u2[,3])
      s1=as.data.frame(prop.table(mytable))
      denom=sum(s1$Freq^q)
      trans=(1/(q-1))*(1-(nem/denom))
    }
  
  feat<-as.matrix(parl)
  feat[fea[1:i]]<-Inf
  idx=which(feat==min(feat[which(feat>0)]))[1]
  fea[1,(i+1)]=idx
}
stopCluster(cl)
 # Feature reduced data with Tsallis entropy
 tsalisdata=data[,fea] 
 return(tsalisdata)
}
