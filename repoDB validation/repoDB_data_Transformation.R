library(data.table)
library(tidyr)
library(purrr)
library(pbapply)
library(reshape2)

repoDB.dat <- fread("Data/repoDB_full.csv")
repoDB.dat <- select(repoDB.dat, c("drug_id", "ind_id", "status"))
repoDB.dat$status = ifelse(repoDB.dat$status == "Approved", 1, 0)
repoDB.dat = tidyr::unite(repoDB.dat, "Ind_Status", c("ind_id", "status"), sep="_")
newDat <- split(repoDB.dat$Ind_Status, repoDB.dat$drug_id)

findStatusofaPair <- function(x,y){
  commonOnes <- intersect(x,y)
  if(length(commonOnes) > 0 & any(endsWith(commonOnes, "_1"))){
    return(1) # TP: If both drugs are approved for the same indication(s)
  }else{
    return(0) # TN: Otherwise
  }
}
repoMat <- pbsapply(newDat, function(x) sapply(newDat, findStatusofaPair, x)) # repoMat is an adjacency matrix (i.e. binary)

repoAdjList <- melt(repoMat) # adjacency matrix to adjacency list (columns: Var1, Var2, value)
repoAdjList <- filter(repoAdjList, value == 1)  # for some reason pipe (%>%) isn't working in my computer
repoAdjList <- filter(repoAdjList, Var1 != Var2)

# ----- get distinct pair only
repoAdjList[1:2] = t(apply(repoAdjList[1:2],1,sort))
repoAdjList = repoAdjList[!duplicated(repoAdjList[1:2]),]
colnames(repoAdjList) <- c("DrugBank.ID1", "DrugBank.ID1", "repoDB_status")

# ----- write into the file
fwrite(repoAdjList, file = "Data/repoDB_drugPair_Transformation.csv")
