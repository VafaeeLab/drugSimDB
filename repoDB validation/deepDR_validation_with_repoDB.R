deepDR_validation_with_repoDB <- function(predDatFile, repoDBFile, threshold=0.01){
  
  require(data.table)
  require(dplyr)
  require(rje)
  require(sets)
  require(pbapply)
  
  # Load prediction data files & get rid of the repeated pairs
  predDat = as.data.frame(fread(file = predDatFile, header = TRUE, sep = ",", quote = "\""))
  predDat = predDat[,c(1,3,5)]
  
  # Load repoDB
  repoDB = fread(file = repoDBFile, header = TRUE, sep = ",", quote = "\"")
  repoDB = repoDB[,c(2,4,6)]
  repoDB[,3] = ifelse(repoDB[,3] == "Approved", 1, 0)
  
  colnames(predDat)[1:2] = colnames(repoDB)[1:2] = c("drug.ID","indication.ID")
  
  # inner join: prediction and gold-standard has the same data-points (based on [drug-disease] composite key)
  joined.table = dplyr::inner_join(predDat,repoDB, by=c("drug.ID","indication.ID"))
  
  # get the score vector
  score = joined.table$`predict score`
  label = joined.table$status
  
  return(list(pred=score,class=label))
  
}
library(ROCit)
library(ggplot2)

val.R.deepDR <- deepDR_validation_with_repoDB(predDatFile = "data/external/deepDR_Table_S10.csv", repoDBFile = "data/repoDB_full.csv")
ROCit_obj = rocit(score=val.R.deepDR$pred, class=val.R.deepDR$class)
# Conclusion: can't run ROC function, as all the predictions were False positives (check with FV )
# TP: Drug-disease that were [approved] as mentioned in the repoDB
# TN: Drug-disease that were [anything but approved] as mentioned in the repoDB
# FP: Drug-disease that were [anything but approved] as mentioned in the repoDB, but predicted by deepDR
# FN: Drug-disease that were [approved] as mentioned in the repoDB, but weren't predicted by deepDR

