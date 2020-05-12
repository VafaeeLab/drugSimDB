library(foreach)
library(doParallel)
library(parallel)
require(pbapply)


setwd("/srv/scratch/z3526914/DrugRepo/")

drug_target <- read.csv("Data/external/uniprot links.csv")
drug_target <- subset(drug_target, Type == "SmallMoleculeDrug")[,c("DrugBank.ID", "UniProt.ID")]
target.similarity <- function(x, y, drug_target){
  tx = as.character(drug_target[which(drug_target[,"DrugBank.ID"]==x),"UniProt.ID"])
  ty = as.character(drug_target[which(drug_target[,"DrugBank.ID"]==y),"UniProt.ID"])
  return(length(intersect(tx, ty))/length(union(tx,ty)))
}



tmp <- unique(subset(drug_target, select = "DrugBank.ID")[,1])



cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("tmp", "target.similarity", "drug_target"), envir = .GlobalEnv)
s <- Sys.time()
sim.matrix<- parallel::parSapply(cl, tmp, function(x) sapply(tmp, target.similarity, x, drug_target))
e <- Sys.time()
print(e-s)



stopCluster(cl)




colnames(sim.matrix) <- tmp
rownames(sim.matrix) <- tmp



# Generate similarity matrix for all SmallMolecul drugs -------------------
drugs.full <- read.csv("Data/external/drug links.csv", row.names = "DrugBank.ID")
drugs.id.full <- rownames(subset(drugs.full, Drug.Type == "SmallMoleculeDrug"))



full.sim.matrix <- matrix(data=NA,nrow=length(drugs.id.full),ncol=length(drugs.id.full))
colnames(full.sim.matrix) <- drugs.id.full
rownames(full.sim.matrix) <- drugs.id.full
full.sim.matrix[rownames(sim.matrix), colnames(sim.matrix)] = sim.matrix



write.csv(full.sim.matrix, paste("Data/ReV_target_Jaccard_similarity.csv", sep =""))