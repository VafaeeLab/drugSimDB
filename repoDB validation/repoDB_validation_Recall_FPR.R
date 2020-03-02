getRepoDB_Status <- function(predDat, repoDB){
  require(data.table)
  require(dplyr)
  require(rje)
  require(sets)
  require(pbapply)
  
  
  # Get the true labels
  getL <- pbapply(predDat, 1, FUN = function(x){
    idx1 = which(repoDB$drug_id %in% x[1])
    idx2 = which(repoDB$drug_id %in% x[2])
    
    # -- if either of drug1 or drug2 is missing in repoDB file,
    if(length(idx1) == 0 | length(idx2) == 0)
      return(NA)
    
    r1 = repoDB[idx1,]
    r2 = repoDB[idx2,]
    
    r1.OR.r2.aprv = (length(r1$status[which(r1$status == "Approved")]) > 0) | (length(r2$status[which(r2$status == "Approved")]) > 0)
    # r2.aprv = r2$status[which(r2$status == "Approved")]
    if(length(intersect(r1$ind_id, r2$ind_id))>0 & r1.OR.r2.aprv) # drug1 and drug2 have similar indication and drug2 is generally an approved drug
      return(1)
    else if (length(intersect(r1$ind_id, r2$ind_id))==0 & r1.OR.r2.aprv) # drug2 is approved, but drug1 and drug2 DON'T have similar indication
      return (0)
    else
      return(0)
    # return(ifelse(length(intersect(r1$ind_id, r2$ind_id))>0 & length(r2.aprv) > 0, 1, 0))
  })
  
  predDat$repoDB_Status = getL
  predDat = predDat[order(predDat$rowMeans, decreasing = T),]
  fwrite(predDat, file="data/new_net_info_V7_pval_RepoDB_Status.csv", na = "NA")
  
  return(predDat)
}

getTopKPrediction <- function(new_predDat, col, from = 20, length.out = 15){
  require(dplyr)
  require(ggplot2)
  new_predDat.2 = new_predDat[-which(is.na(new_predDat$repoDB_Status)),]
  colIndx = grep(col,colnames(new_predDat))
  new_predDat.2 = new_predDat.2[,c("ID1","ID2",col, "repoDB_Status")]
  new_predDat.2 = new_predDat.2[order(new_predDat.2[,3], decreasing = T),]
  if(length(which(is.na(new_predDat.2[,3]))) != 0){
    new_predDat.2 = new_predDat.2[-which(is.na(new_predDat.2[,3])),]
  }
  
  retDat = NULL
  breaks = round(exp(seq(log(from), log(nrow(new_predDat)), length.out = length.out)))
  # breaks = round(exp(seq(log(from), log(10000), length.out = length.out)))
  
  for(k in breaks){
    # print(k)
    new_predDat.filt = new_predDat.2[1:k, ]
    nTP = nrow(new_predDat.filt[which(new_predDat.filt$repoDB_Status == 1),])
    nFP = nrow(new_predDat.filt[which(new_predDat.filt$repoDB_Status == 0),])
    
    recall_at_k = nTP / k
    fpr_at_k = nFP / k
    # print(recall_at_k)
    # precision_at_k = 1 - recall_at_k
    retDat = rbind(retDat, cbind(k, recall_at_k, fpr_at_k))
  }
  retDat = as.data.frame(retDat)
  colnames(retDat) = c("k", paste0(col,".Recall"), paste0(col,".FPR"))
  # p = ggplot(data=retDat, 
  #            aes(x=k, y=recall_at_k, group=1)) +
  #     geom_line(color = "red", linetype="dashed") +
  #     geom_point()
  # return(p)
  return(retDat)
}


library(ROCit)
library(ggplot2)
library(dplyr)

predDatFile = "data/new_net_info_V7_pval.csv"
repoDBFile = "data/repoDB_full.csv"
predDat = as.data.frame(fread(file = predDatFile, header = TRUE, sep = ",", quote = "\""))
predDat[1:2] = t(apply(predDat[1:2],1,sort))
predDat = predDat[!duplicated(predDat[1:2]),]

repoDB = fread(file = repoDBFile, header = TRUE, sep = ",", quote = "\"")

new_predDat <- getRepoDB_Status(predDat, repoDB)

col = "Chem_similarity"
length.out = 50

retDat <- getTopKPrediction(new_predDat,col = col, from = 20, length.out = length.out)
p1 = ggplot(data=retDat, aes(x=k, y=paste0(col,".Recall"), group=1)) +
  geom_line(color = "red", linetype="dashed") + geom_point()
ggsave(paste0("topk_Recall@k_",col,".pdf"),
       p1, device = "pdf", path = "data/")
p2 = ggplot(data=retDat, aes(x=k, y=paste0(col,".FPR"), group=1)) +
  geom_line(color = "red", linetype="dashed") + geom_point()
ggsave(paste0("topk_FPR@k_",col,".pdf"),
       p2, device = "pdf", path = "data/")

col = "Target_similarity"
newRetDat <- getTopKPrediction(new_predDat,col = col, from = 20, length.out = length.out)
p1 = ggplot(data=newRetDat, aes(x=k, y=paste0(col,".Recall"), group=1)) +
  geom_line(color = "red", linetype="dashed") + geom_point()
ggsave(paste0("topk_Recall@k_",col,".pdf"),
       p1, device = "pdf", path = "data/")
p2 = ggplot(data=newRetDat, aes(x=k, y=paste0(col,".FPR"), group=1)) +
  geom_line(color = "red", linetype="dashed") + geom_point()
ggsave(paste0("topk_FPR@k_",col,".pdf"),
       p2, device = "pdf", path = "data/")

retDat = dplyr::inner_join(retDat,newRetDat, by="k")

col = "Pathway_similarity"
newRetDat <- getTopKPrediction(new_predDat,col = col, from = 20, length.out = length.out)
p1 = ggplot(data=newRetDat, aes(x=k, y=paste0(col,".Recall"), group=1)) +
  geom_line(color = "red", linetype="dashed") + geom_point()
ggsave(paste0("topk_Recall@k_",col,".pdf"),
       p1, device = "pdf", path = "data/")
p2 = ggplot(data=newRetDat, aes(x=k, y=paste0(col,".FPR"), group=1)) +
  geom_line(color = "red", linetype="dashed") + geom_point()
ggsave(paste0("topk_FPR@k_",col,".pdf"),
       p2, device = "pdf", path = "data/")
retDat = dplyr::inner_join(retDat,newRetDat, by="k")

col = "GO_CC_Similarity"
newRetDat <- getTopKPrediction(new_predDat,col = col, from = 20, length.out = length.out)
p1 = ggplot(data=newRetDat, aes(x=k, y=paste0(col,".Recall"), group=1)) +
  geom_line(color = "red", linetype="dashed") + geom_point()
ggsave(paste0("topk_Recall@k_",col,".pdf"),
       p1, device = "pdf", path = "data/")
p2 = ggplot(data=newRetDat, aes(x=k, y=paste0(col,".FPR"), group=1)) +
  geom_line(color = "red", linetype="dashed") + geom_point()
ggsave(paste0("topk_FPR@k_",col,".pdf"),
       p2, device = "pdf", path = "data/")
retDat = dplyr::inner_join(retDat,newRetDat, by="k")

col = "GO_MF_Similarity"
newRetDat <- getTopKPrediction(new_predDat,col = col, from = 20, length.out = length.out)
p1 = ggplot(data=newRetDat, aes(x=k, y=paste0(col,".Recall"), group=1)) +
  geom_line(color = "red", linetype="dashed") + geom_point()
ggsave(paste0("topk_Recall@k_",col,".pdf"),
       p1, device = "pdf", path = "data/")
p2 = ggplot(data=newRetDat, aes(x=k, y=paste0(col,".FPR"), group=1)) +
  geom_line(color = "red", linetype="dashed") + geom_point()
ggsave(paste0("topk_FPR@k_",col,".pdf"),
       p2, device = "pdf", path = "data/")
retDat = dplyr::inner_join(retDat,newRetDat, by="k")

col = "GO_BP_Similarity"
newRetDat <- getTopKPrediction(new_predDat,col = col, from = 20, length.out = length.out)
p1 = ggplot(data=newRetDat, aes(x=k, y=paste0(col,".Recall"), group=1)) +
  geom_line(color = "red", linetype="dashed") + geom_point()
ggsave(paste0("topk_Recall@k_",col,".pdf"),
       p1, device = "pdf", path = "data/")
p2 = ggplot(data=newRetDat, aes(x=k, y=paste0(col,".FPR"), group=1)) +
  geom_line(color = "red", linetype="dashed") + geom_point()
ggsave(paste0("topk_FPR@k_",col,".pdf"),
       p2, device = "pdf", path = "data/")
retDat = dplyr::inner_join(retDat,newRetDat, by="k")

col = "rowMeans"
newRetDat <- getTopKPrediction(new_predDat,col = col, from = 20, length.out = length.out)
p1 = ggplot(data=newRetDat, aes(x=k, y=paste0(col,".Recall"), group=1)) +
  geom_line(color = "red", linetype="dashed") + geom_point()
ggsave(paste0("topk_Recall@k_",col,".pdf"),
       p1, device = "pdf", path = "data/")
p2 = ggplot(data=newRetDat, aes(x=k, y=paste0(col,".FPR"), group=1)) +
  geom_line(color = "red", linetype="dashed") + geom_point()
ggsave(paste0("topk_FPR@k_",col,".pdf"),
       p2, device = "pdf", path = "data/")
retDat = dplyr::inner_join(retDat,newRetDat, by="k")

# all together
p1 <- ggplot(retDat, aes(x=k)) + 
  geom_line(aes(y = Chem_similarity.Recall), color="steelblue", linetype="twodash", ) + 
  geom_line(aes(y = Target_similarity.Recall), color="purple", linetype="twodash") + 
  geom_line(aes(y = Pathway_similarity.Recall), color="cyan", linetype="twodash") + 
  geom_line(aes(y = GO_CC_Similarity.Recall), color="green", linetype="twodash") + 
  geom_line(aes(y = GO_MF_Similarity.Recall), color="black", linetype="twodash") + 
  geom_line(aes(y = GO_BP_Similarity.Recall), color="yellow", linetype="twodash") + 
  geom_line(aes(y = rowMeans.Recall), color="red", linetype="twodash") +
  labs(y="Recall@topK", x="topK")
ggsave(paste0("topk_all_Recall@K_",col,".pdf"),
       p1, device = "pdf", path = "data/")

p2 <- ggplot(retDat, aes(x=k)) + 
  geom_line(aes(y = Chem_similarity.FPR), color="steelblue", linetype="twodash", ) + 
  geom_line(aes(y = Target_similarity.FPR), color="purple", linetype="twodash") + 
  geom_line(aes(y = Pathway_similarity.FPR), color="cyan", linetype="twodash") + 
  geom_line(aes(y = GO_CC_Similarity.FPR), color="green", linetype="twodash") + 
  geom_line(aes(y = GO_MF_Similarity.FPR), color="black", linetype="twodash") + 
  geom_line(aes(y = GO_BP_Similarity.FPR), color="yellow", linetype="twodash") + 
  geom_line(aes(y = rowMeans.FPR), color="red", linetype="twodash") +
  labs(y="FPR@topK", x="topK")
ggsave(paste0("topk_all_FPR@K_",col,".pdf"),
       p2, device = "pdf", path = "data/")