install.packages('digest')
install.packages('rlang')
install.packages('biomaRt')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")
install.packages("remotes")
remotes::install_github("gongx030/scDatasets")


library(SingleCellExperiment)
library(edgeR)
library(scDatasets)
library(biomaRt)
library('Linnorm')
normalized_data =  function (raw_Data, min_Reads = 5, min_Cell = 0.1, min_Gene = 1000,
                             log = T)
{
  
  #cell Filtering
  MIN = min(raw_Data)
  good_Cells = apply(raw_Data, 2, function(x) sum(x > MIN)) >= min_Gene
  temp_Data = raw_Data[, good_Cells]
  
  #Gene filtering
  C = floor(dim(temp_Data)[2] * min_Cell)
  exprs_Genes = apply(temp_Data, 1, function(x) sum(x > min_Reads)) >= C
  temp_Data = temp_Data[exprs_Genes, ]
  
  cat(paste("Remaining genes:", dim(temp_Data)[1], "\n", sep = ""))
  cat(paste("Remaining cells:", dim(temp_Data)[2], sep = ""))
  norm_Data = Linnorm.Norm(temp_Data)
  if (log == T) {
    norm_log_Data = log2(1 + norm_Data)
  }
  else {
    norm_log_Data = norm_Data
  }
  return(norm_log_Data)
}


##################Biase_data###########################################
#Reading RDS files for rds file (darmanis, yan) and download pollen data from scDatasets
Biase_data<- readRDS("yan.rds")
#Baise_data=scDatasets::pollen
data <- assay(Biase_data) 

annotation <- Biase_data[[1]] #already factor type class
colnames(data) <- annotation
Normalized_Biase_data = normalized_data(datan)


#saving file in csv format, cells in row, genes in coloumn
write.table(t(Normalized_Biase_data),file="preprocessdata.csv",sep=",")
