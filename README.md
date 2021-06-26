# sc-REnF

An Entropy Based Feature Selection- Application on Single cell RNA Sequence Data

## Pre-requisites

> R version  4.0.0 or higher

> Python 3.7

## Install
library("devtools")

install_github("Snehalikalall/sc-REnF")

Check the installation:

library(scREnF)

## Load required packages

R packages

     library(SingleCellExperiment)
     library(edgeR)
     library(MASS)
     library(foreach)
     library(doParallel)
     library(biomaRt)
     library(Linnorm)
     library(splatter)

Python Packages: 
 
    pip install scanpy
    pip install leidenalg


## Usage of the R functions

Preprocess raw data using normalized_data.R function

    Biase_data<- readRDS("yan.rds")
    data <- assay(Biase_data) 
    annotation <- Biase_data[[1]] #already factor type class
    colnames(data) <- annotation
    preprocessed_data = normalized_data(data, 5, 0.1, 1000, T)


Use Renyifeature.R to select features using Renyi entropy, Tsallisfeature.R to select features using Tsallis entropy

    # load the preprocess data
    data=t(as.matrix(read.csv("yan_process.csv",header=FALSE)))
    cell<-as.matrix(read.csv("celltype.csv",header=FALSE))
    n <- nrow(data)
    col<-ncol(data)
    count=ncol(data)
    #nf: Number of feature to be selected, default is 500; P: Number of cores, default is 30;
    # q-value = 0.3 for tsallis, q-value=0.7 for Renyi

    # Renyi entropy based Feature Selection, the function returns data with selected features
    Feadata=Renyifeature(data,cell,p,q,nf)
    
    # Renyi entropy based Feature Selection, the function returns the data with selected features
    Feadata=Tsallisfeature(data,cell,p,q,nf)



# A dry run on CBMC data 

[A demo run of sc-REnF for CBMC data](https://snehalikalall.github.io/Introduction-to-sc-REnF/)


