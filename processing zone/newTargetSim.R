library("pbapply")

# ### Option list
args = commandArgs(trailingOnly=FALSE)

### Make Cluster nodes for parallelizing
dataDir = args[4]
protSimFile = args[5]
drugTargetFile = args[6]
startRow = strtoi(args[7])
offset = strtoi(args[8])
outFile = args[9]

print(args)


prot.sim.matrix <- read.csv(paste0(dataDir, protSimFile))

getIDs <- function(str){
  return(gsub("drugbank_target.","", str))
}
protIds <- sapply(as(colnames(prot.sim.matrix), "list"), getIDs)
colnames(prot.sim.matrix) <- protIds
rownames(prot.sim.matrix) <- protIds
prot.sim.matrix = prot.sim.matrix/100

# Read drug and their targets -------------------------------------------------------
# drug_target <- read.csv("/short/yr31/aa7970/azData/DD/Data/uniprot links.csv")
drug_target <- read.csv(paste0(dataDir, drugTargetFile))
drug_target <- subset(drug_target, Type == "SmallMoleculeDrug")[,c("DrugBank.ID", "UniProt.ID")]

# compute target-based drug similarity -------------------------------------------------------

# this function compute similarity between drug x and y based on their targets
target.similarity <- function(x, y, drug_target, prot.sim.matrix){
  tx = as.character(drug_target[which(drug_target[,"DrugBank.ID"]==x),"UniProt.ID"])
  ty = as.character(drug_target[which(drug_target[,"DrugBank.ID"]==y),"UniProt.ID"])
  return(mean(prot.sim.matrix[tx,ty]))
}

tmp <- unique(subset(drug_target, select = "DrugBank.ID")[,1])

endRow = startRow + offset -1
endRow = if (endRow > length(tmp)) length(tmp) else endRow


st <- Sys.time()
sim.matrix<- pbsapply(tmp, function(x) sapply(tmp[startRow:endRow], target.similarity, x, drug_target, prot.sim.matrix))
write.csv(sim.matrix, paste0(dataDir, outFile,"_",startRow, "_",endRow,".csv"), quote = F, row.names = F)
e<- Sys.time()
print(e-st)

# write.csv(sim.matrix, paste("/short/yr31/aa7970/azData/DD/Data/Target_similarity.csv", sep =""))

