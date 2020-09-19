library('MASS')
library(foreach)
library(doParallel)
library('Rfast')

# load the preprocess data
data=as.matrix(read.csv("preprocessdata.csv"))
#Cell type vector
cell<-rownames(data)  
n <- nrow(data)
col<-ncol(data)
count=ncol(data)            
#Number of feature to be selected,nf=50 default
# Number of cores, p= 20 default
p=20
# q-value =0.3 default for tsallis, q-value=0.7 for default Renyi
# Renyi entropy based Feature Selection
Renyifeaturedata=function(data,cell,p,q,nf){
  jointprob_2 <-function(u) {
    ncol=length(unique(u[,1]))
    nrow=length(unique(u[,2]))
    probtable=matrix(0,nrow,ncol)
    a=unique(u[,1])
    for(i in 1:length(a))
    {
      index=which(u[,1]==a[i])
      pmat1=as.matrix(u[index,2])
      s1=aggregate(data.frame(count = pmat1), list(value = pmat1), length)
      lens1=length(s1$count)
      probtable[1:lens1,i]=(s1$count)/nrow(data)
    }  
    probtable
  }
  jointprob_3 <-function(u) {
    ncol=length(unique(u[,1]))
    nrow=length(unique(u[,2]))*length(unique(u[,3]))
    probtable=matrix(0,nrow,ncol)
    a=unique(u[,1])
    for(i in 1:length(a))
    {
      index=which(u[,1]==a[i])
      pmat=as.matrix(u[index,2:3])
      b=unique(pmat[,1])
      count=1
      for(j in 1: length(b))
      {
        index1=which(pmat[,1]==b[j])
        pmat1=pmat[index1,2]
        s1=aggregate(data.frame(count = pmat1), list(value = pmat1), length)
        lens1=length(s1$count)
        probtable[(count:(count+lens1-1)),i]=(s1$count)/nrow(data)
        count=count+lens1
      }
    }
    probtable
  }
  
  
  cl <- makeCluster(p) 
  registerDoParallel(cl)
  classc=cell
  fea<- matrix(0, nrow=1,ncol =nf)
  parl<-foreach(j=1:count, .combine=c,.packages= c("Rfast")) %dopar%
  {
    u<-as.matrix(cbind(classc,data[,j]))
    tmat=jointprob_2(u)
    tmat=tmat^q
    transvecfea=(q/(1-q))*(log(sum((rowsums(tmat))^(1/q))))
  }
  
  mimc1<-as.matrix(parl)
  fea[1,1]=which.max(mimc1)
  idx=which.max(mimc1)
 
  for(i in 1:(nf-1))
  {  
    parl<-foreach(j=1:count, .combine=c,.packages= c("Rfast")) %dopar%
    {
      u2<-as.matrix(cbind(classc,data[,idx],data[,j]))
      tmat=jointprob_3(u2)
      tmat=tmat^q
      transvecfea=(q/(1-q))*(log(sum((rowsums(tmat))^(1/q))))
    }
    feat<-as.matrix(parl)
    feat[fea[1:i]]<-Inf
    idx=which(feat==min(feat[which(feat>0)]))[1]
    fea[1,(i+1)]=idx
  }
  stopCluster(cl)
  Renyidata=data[,fea] # Feature reduced Data with renyi entropy 
}

# Tsallis Entropy based feature selection
Tsallisfeaturedata=function(data,cell,p,q,nf){
  jointprob_2 <-function(u) {
    ncol=length(unique(u[,1]))
    nrow=length(unique(u[,2]))
    probtable=matrix(0,nrow,ncol)
    a=unique(u[,1])
    for(i in 1:length(a))
    {
      index=which(u[,1]==a[i])
      pmat1=as.matrix(u[index,2])
      s1=aggregate(data.frame(count = pmat1), list(value = pmat1), length)
      lens1=length(s1$count)
      probtable[1:lens1,i]=(s1$count)/nrow(data)
    }  
    probtable
  }
  jointprob_3 <-function(u) {
    ncol=length(unique(u[,1]))
    nrow=length(unique(u[,2]))*length(unique(u[,3]))
    probtable=matrix(0,nrow,ncol)
    a=unique(u[,1])
    for(i in 1:length(a))
    {
      index=which(u[,1]==a[i])
      pmat=as.matrix(u[index,2:3])
      b=unique(pmat[,1])
      count=1
      for(j in 1: length(b))
      {
        index1=which(pmat[,1]==b[j])
        pmat1=pmat[index1,2]
        s1=aggregate(data.frame(count = pmat1), list(value = pmat1), length)
        lens1=length(s1$count)
        probtable[(count:(count+lens1-1)),i]=(s1$count)/nrow(data)
        count=count+lens1
      }
    }
    probtable
  }
  
  cl <- makeCluster(p)
  registerDoParallel(cl)
  classc=cell
  fea<- matrix(0, nrow=1,ncol =nf)
  parl<-foreach(j=1:count, .combine=c,.packages= c("Rfast")) %dopar%
  {
    u<-as.matrix(cbind(classc,data[,j]))
    tmat=jointprob_2(u)
    tmatq=tmat^q
    nem=sum(tmatq)
    s1=aggregate(data.frame(count = data[,j]), list(value = data[,j]), length)
    lens1=length(s1$count)
    probnew=(s1$count)/nrow(data)
    denom=sum(probnew^q)
    trans=(1/(q-1))*(1-(nem/denom))
  }
  
  mimc1<-as.matrix(parl)
  fea[1,1]=which.max(mimc1)
  idx=which.max(mimc1)
  #fea[1,1]<-which(mimc1==max(mimc1[which(mimc1>0)]))[1]
  #idx=which(mimc1==min(mimc1[which(mimc1>0)]))[1]
  
  for(i in 1:(nf-1))
  {  
    parl<-foreach(j=1:count, .combine=c,.packages= c("Rfast")) %dopar%
    {
      u2<-as.matrix(cbind(classc,data[,idx],data[,j]))
      tmat=jointprob_3(u2)
      tmatq=tmat^q
      nem=sum(tmatq)
      u3<-as.matrix(cbind(data[,idx],data[,j]))
      probnew=jointprob_2(u3)
      denom=sum(probnew^q)
      trans=(1/(q-1))*(1-(nem/denom))
    }
    feat<-as.matrix(parl)
    feat[fea[1:i]]<-Inf
    idx=which(feat==min(feat[which(feat>0)]))[1]
    fea[1,(i+1)]=idx
  }
  tsalisdata=data[,fea] # Feature reduced data with Tsallis data
  
}

#Run the code
Feadata=Tsallisfeaturedata(data,cell,p,q,nf)


#Saving feature selected data for further analysis
write.table(Feadata,file="DatafeaSelected.csv",sep=",")

