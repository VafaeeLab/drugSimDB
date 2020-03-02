
# KEGG pathways associated with each drug: Uncomment if update is needed -------------------------------------------
# source("./update_KEGGpathways_for_drugs.R")
# update_KEGGpathways_for_drugs()

# Compute pathway similarities ------------------------------------------------------------------------
setwd("/srv/scratch/z3526914/DrugRepo/")
library(pbapply)
library(data.table)
library(dplyr)
library("org.Hs.eg.db")

drug_path <- read.csv("Data/durg_KEGGpathways.csv")
drug_path[,2] <- gsub("ko","",drug_path$Kegg.ID,fixed = T)

drug2Path <- split(drug_path[,2], drug_path[,1])
message("Loading Drug-KEGG mapping file: [DONE]")

# Get KEGG genes  
genesKegg <- as.list(org.Hs.egPATH)
message("Loading KEGG genes: [DONE]")

path.similarity <- function(d1,d2){
  # path.set1 <- unlist(d1[],use.names = F); path.set2 <- unlist(d2[],use.names = F);
  path.set1 <- drug2Path[[d1]]; path.set2 <- drug2Path[[d2]]; 
  val <- lapply(path.set1,FUN = function(x){
    list1 <- lapply(path.set2, FUN = function(y){
      return(pathSim(x,y,genesKegg))
    })
    return(max(unlist(list1), na.rm = T))
  }) 
  return(max(unlist(val), na.rm = T))
}

message("Computing non sim.matrix")
s <- Sys.time()

# ---------- non-parallel execution
# sim.matrix <- pbsapply(seq_len(length(drug2Path)), function(x) sapply(seq_len(length(drug2Path)), path.similarity, x))

#----------- parallel excution
library(parallel)
numCores <- detectCores()
cl <- makeCluster(numCores)
clusterEvalQ(cl, {
  library(BioCor)
})
clusterExport(cl=cl, varlist=c("genesKegg","drug2Path", "path.similarity"), envir = .GlobalEnv)

sim.matrix <- parallel::parSapply(cl, seq_len(length(drug2Path)), function(x) sapply(seq_len(length(drug2Path)), path.similarity, x))
stopCluster(cl)

e <- Sys.time()
message("[DONE]")
print(e-s)

# Generate similarity matrix for all SmallMolecul drugs -------------------
drugs.full <- read.csv("Data/external/drug links.csv", row.names = "DrugBank.ID")
drugs.id.full <- rownames(subset(drugs.full, Drug.Type == "SmallMoleculeDrug"))
message("Loading drug links file: [DONE]")

full.sim.matrix <- matrix(data=NA,nrow=length(drugs.id.full),ncol=length(drugs.id.full))
colnames(full.sim.matrix) <- drugs.id.full
rownames(full.sim.matrix) <- drugs.id.full
full.sim.matrix[rownames(sim.matrix), colnames(sim.matrix)] = sim.matrix
message("Generating full sim matrix: [DONE]")

# Save the similarity matrix ----------------------------------------------
fwrite(as.data.frame(full.sim.matrix), paste("Data/","path_similarity_v2",".csv", sep =""), na = "NA", row.names = T)




