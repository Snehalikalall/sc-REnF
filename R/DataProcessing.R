normalized_data <- function (raw_Data, min_Reads = 5, min_Cell = 0.1, min_Gene = 1000,
                             log = T)
{
 library(SingleCellExperiment)
  library(edgeR)
  library('Linnorm')

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
