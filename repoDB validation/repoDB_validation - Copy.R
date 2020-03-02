repoDB_validation <- function(predDatFile, repoDBFile, threshold){
  require(data.table)
  require(dplyr)
  require(rje)
  require(sets)
  require(pbapply)

  # Load data files
  predDat = fread(file = predDatFile, header = TRUE, sep = ",", quote = "\"")
  repoDB = fread(file = repoDBFile, header = TRUE, sep = ",", quote = "\"")
  
  # Filter predicted Drug-pairs for which drug data available in repoDB
  predDat = predDat[which((predDat$ID1 %in% repoDB$drug_id) & (predDat$ID2 %in% repoDB$drug_id)),]
  predDat = predDat[which(predDat$p_value < 0.01),]
  
  # # Get TP and TN info from repoDB
  # generate Predicted labels
  # P_vector = ifelse(predDat$adjP_value < threshold, 1, 0)
  # OR, use "raw" p-values
  P_vector = 1 - predDat$adjP_value
  getL <- pbapply(predDat, 1, FUN = function(x){
    r1 = repoDB[which(repoDB$drug_id %in% x[1]),]
    r2 = repoDB[which(repoDB$drug_id %in% x[2]),]
    r2.aprv = r2$ind_id[which(r2$status == "Approved")]
    if(length(intersect(r1$ind_id, r2$ind_id))>0 & length(r2.aprv) > 0)
      return(1)
    else if (length(intersect(r1$ind_id, r2$ind_id))==0 & length(r2.aprv) > 0)
      return (NA)
    else
      return(0)
    # return(ifelse(length(intersect(r1$ind_id, r2$ind_id))>0 & length(r2.aprv) > 0, 1, 0))
  })
  
  return(list(pred=P_vector,class=getL))
}

library(ROCit)
library(ggplot2)
# for full Chemical similarity
# val.R <- repoDB_validation(predDatFile = "data/",repoDBFile = "data/repoDB_full.csv", threshold = 0.01)
# ROCit_obj = rocit(score=val.R$pred, class=val.R$class)
# plot(ROCit_obj)


# for full data
val.R <- repoDB_validation(predDatFile = "data/new_net_info_V4_pval.csv",repoDBFile = "data/repoDB_full.csv", threshold = 0.01)
newPred = val.R$pred[-which(is.na(val.R$class))]
newClass = val.R$class[-which(is.na(val.R$class))]
ROCit_obj = rocit(score=newPred, class=newClass)
plot(ROCit_obj)
