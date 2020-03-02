library(gtools)
ontoTag = "BP"
datadir <- paste0("Data/GO_Sim/",ontoTag, "/")
files <- mixedsort(list.files(path = datadir, pattern="*.csv", recursive=FALSE, full.names=TRUE), decreasing = F)

mat <- NULL
for(i in 1:length(files)){
  # for(i in 1:2){
  temp <- read.table(file = files[i], header = T, sep = ",")
  mat = rbind(mat,temp)
}
# fix the column header
# Drugs names that were used for this similarity measurement
drugs.full <- read.csv("Data/external/drug links.csv", row.names = "DrugBank.ID")
drugs.id.full <- rownames(subset(drugs.full, Drug.Type == "SmallMoleculeDrug"))

colnames(mat) <- drugs.id.full
rownames(mat) <- drugs.id.full
library(data.table)
fwrite(mat, paste0(datadir, "GO_Sim_", ontoTag, "_combined.csv"), row.names = T, na = "NA")
# rm(mat)
