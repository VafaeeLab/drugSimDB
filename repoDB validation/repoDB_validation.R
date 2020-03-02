repoDB_validation <- function(predDatFile, repoDBFile, threshold){
  require(data.table)
  require(dplyr)
  require(rje)
  require(sets)
  require(pbapply)

  # Load prediction data files & get rid of the repeated pairs
  predDat = as.data.frame(fread(file = predDatFile, header = TRUE, sep = ",", quote = "\""))
  # predDat[1:2] = t(apply(predDat[1:2],1,sort))
  # predDat = predDat[!duplicated(predDat[1:2]),]
  # save 
  fwrite(predDat, file="data/repoDB_PredDat_V7.csv", sep=",", row.names = F, quote = T)
  
  repoDB = fread(file = repoDBFile, header = TRUE, sep = ",", quote = "\"")
  
  # Filter predicted Drug-pairs for which drug data available in repoDB
  predDat = predDat[which((predDat$ID1 %in% repoDB$drug_id) & (predDat$ID2 %in% repoDB$drug_id)),]
  # predDat = predDat[which(predDat$p_value < 0.01),]
  
  # # Get TP and TN info from repoDB
  # generate Predicted labels
  # P_vector = ifelse(predDat$adjP_value < threshold, 1, 0)
  # OR, use "raw" p-values
  
  # P_vector = 1 - predDat$adjP_value
  P_vector = predDat$rowMeans # this isn't appropriate, because, its just a combined score, doesn't indicate probability
  
  # Get the true labels
  getL <- pbapply(predDat, 1, FUN = function(x){
    r1 = repoDB[which(repoDB$drug_id %in% x[1]),]
    r2 = repoDB[which(repoDB$drug_id %in% x[2]),]
    r1.OR.r2.aprv = (length(r1$status[which(r1$status == "Approved")]) > 0) | (length(r2$status[which(r2$status == "Approved")]) > 0)
    # r2.aprv = r2$status[which(r2$status == "Approved")]
    if(length(intersect(r1$ind_id, r2$ind_id))>0 & r1.OR.r2.aprv) # drug1 and drug2 have similar indication and drug2 is generally an approved drug
      return(1)
    else if (length(intersect(r1$ind_id, r2$ind_id))==0 & r1.OR.r2.aprv) # drug2 is approved, but drug1 and drug2 DON'T have similar indication
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
val.R <- repoDB_validation(predDatFile = "Data/new_net_info_V7_pval.csv",repoDBFile = "Data/repoDB_full.csv", threshold = 0.01)
newPred = val.R$pred[-which(is.na(val.R$class))]
newClass = val.R$class[-which(is.na(val.R$class))]

# ------------ Simple ROC 
ROCit_obj = rocit(score=newPred, class=newClass)
ROCit_obj$AUC

# ------------ ROC + confindence interval
score = newPred
class = newClass
rocit_emp <- rocit(score = score, 
                   class = class, 
                   method = "emp")
rocit_bin <- rocit(score = score, 
                   class = class, 
                   method = "bin")
# --------------------------
ciROC_emp90 <- ciROC(rocit_emp, 
                     level = 0.9)
set.seed(200)
ciROC_bin90 <- ciROC(rocit_bin, 
                     level = 0.9, nboot = 200)
plot(ciROC_emp90, col = 1, 
     legend = FALSE)
lines(ciROC_bin90$TPR~ciROC_bin90$FPR, 
      col = 2, lwd = 2)
lines(ciROC_bin90$LowerTPR~ciROC_bin90$FPR, 
      col = 2, lty = 2)
lines(ciROC_bin90$UpperTPR~ciROC_bin90$FPR, 
      col = 2, lty = 2)
legend("bottomright", c("Empirical ROC",
                        "Binormal ROC",
                        "90% CI (Empirical)", 
                        "90% CI (Binormal)"),
       lty = c(1,1,2,2), col = 
         c(1,2,1,2), lwd = c(2,2,1,1))

# ----------------- KS plot
# KS plot shows the cumulative density functions F(c) and G(c) in the positive 
# and negative populations. If the positive population have higher value, then 
# negative curve (F(c)) ramps up quickly. The KS statistic is the maximum difference 
# of F(c) and G(c). (Source: https://cran.r-project.org/web/packages/ROCit/vignettes/my-vignette.html)
kplot <- ksplot(ROCit_obj)
message("KS Stat (empirical) : ", kplot$`KS stat`)
message("KS Stat (empirical) cutoff : ", kplot$`KS Cutoff`)

# ---------------- Gain table
# For description: https://cran.r-project.org/web/packages/ROCit/vignettes/my-vignette.html)
rocit_emp <- rocit(score = score, 
                   class = class, 
                   negref = "FP")
gtable_custom <- gainstable(rocit_emp, 
                            breaks = seq(1,100,15))
plot(gtable15, type = 1)

# ----------------- Precision-vs-Recall curve
# ACC: Overall accuracy of classification.
# MIS: Misclassification rate.
# SENS: Sensitivity.
# SPEC: Specificity.
# PREC: Precision.
# REC: Recall. Same as sensitivity.
# PPV: Positive predictive value.
# NPV: Positive predictive value.
# TPR: True positive rate.
# FPR: False positive rate.
# TNR: True negative rate.
# FNR: False negative rate.
# pDLR: Positive diagnostic likelihood ratio.
# nDLR: Negative diagnostic likelihood ratio.
# FSCR: F-score,.

measure <- measureit(score = score, class = class,
                     measure = c("PREC", "REC", "FSCR", "ACC"))
names(measure)
plot(measure$PREC~measure$REC, type = "l")
measure$FSCR
plot(measure$ACC~measure$Cutoff, type = "l")
