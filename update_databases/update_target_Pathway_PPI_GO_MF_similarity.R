require(igraph)
require(GOSemSim)
require(pbapply)
require(clusterProfiler)
require(org.Hs.eg.db)
require(stringi)
require(AnnotationDbi)

setwd("/srv/scratch/z3526914/DrugRepo/")
# 1. Download GO terms from EnrichR library (https://amp.pharm.mssm.edu/Enrichr/#stats)
# MF
pop.filepath = "Data/external/Enrichr/GO_Molecular_Function_2018.txt"
num.col = max(count.fields(pop.filepath, sep = "\t", blank.lines.skip = TRUE), na.rm = TRUE)
GO_MF = read.table(pop.filepath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE, col.names = 1:num.col, blank.lines.skip = TRUE, quote = "")
# filter specific-general terms
nGS.terms.MF = rowSums(GO_MF != "", na.rm = T)-1
newMF = GO_MF[-which(rowSums(GO_MF != "", na.rm = T)-1 <= 15),]
newMF = newMF[-which(rowSums(newMF != "", na.rm = T)-1 >= 100),]
newMF = dplyr::select(newMF, -2) # throw away the the column with all NA (2nd-Column)

# 2. 
all.GO.terms = newMF
# Parse the GO accession numbers for each Terms
row.names(all.GO.terms) = gsub(")","",lapply(lapply(stri_split_fixed(stri_reverse(all.GO.terms[,1]),"(",n = 2), FUN = stri_reverse), "[[",1))
nTerms = nrow(all.GO.terms)
pop.genes = unique(c(as.matrix(all.GO.terms[,-1])))
N = length(pop.genes)

# 
# Download Human PPI network from i2d database
i2d.db = read.delim("Data/external/i2d.2_9.Public.HUMAN.tab", sep="\t", header = T)
# Load the PPI as an iGraph object and remove cycles and loops
ppi.net.all = igraph::simplify(graph.data.frame(i2d.db[,c("SwissProt1","SwissProt2")], directed = F))
rm(i2d.db)

# The drug-pairs that matters
drugRepo.db = read.delim("Data/new_net_info_V3_pval.csv", sep=",", header = T)
# The target proteins of each Drugs
drugTarget.db = read.delim("Data/external/uniprot links.csv", sep=",", header = T)
# only focus on the Small Molecule drugs
drugTarget.db = drugTarget.db[which(drugTarget.db$Type == "SmallMoleculeDrug"),c(1,4)]
# split the db to get the target grouped for each drugs: it becomes a list now
drugTarget.list = split(as.character(drugTarget.db[,2]), as.character(drugTarget.db[,1]))
rm(drugTarget.db)


drugs.full <- read.csv("Data/external/drug links.csv", row.names = "DrugBank.ID")
drugs.id.full <- rownames(subset(drugs.full, Drug.Type == "SmallMoleculeDrug"))

allTargetList = lapply(drugs.id.full,FUN = function(x){
  indx = which(names(drugTarget.list) == x)
  if(length(indx)==0)   # if the drug doesn't have any target info
    return("")
  else
    return(drugTarget.list[[indx]])
})
# make named list (i.e. a dictionary where key: drugId, and value: targetlist)
names(allTargetList) = drugs.id.full

# Get PPI neighbours of a Protein
getNeighbours <- function(aNode){
  #if the node doesn't exist in the net, return null
  if(sum(which(V(ppi.net.all)$name == aNode)) == 0)
    return(aNode)
  nb <- as.character(ego(ppi.net.all, order=1, nodes = aNode, mode = "all", mindist = 0)[[1]]$name)
  if(length(nb) == 1){# indicates no PPI parterns exists, except itself (since mindist=0)
    return(aNode)
  }else{  # return the PPI partners except itself
    # return(setdiff(nb,aNode))
    return (nb)
  }
}
getGOTerms <- function(i){
  info.gene.overrep = data.frame(matrix(Inf, nrow = nrow(all.GO.terms), ncol = 3))
  # get all the targets associated with the aDrug
  # t.prots = aDrug[]
  t.prots = allTargetList[[i]]
  # find a Set of PPI neighbours of all the target proteins
  t.prots.neigh = unique(unlist(lapply(unlist(t.prots), 
                                       FUN = function(x){return(getNeighbours(x))})))
  # Get Gene Symbols of each UniProt Proteins
  t.gns.neigh <- tryCatch(expr = intraIDMapper(t.prots.neigh, species = "HOMSA", 
                                               srcIDType = "UNIPROT", destIDType = "SYMBOL"),
                          error = function(e) return(""))
  t.gns.neigh <- t.gns.neigh %>% base::unlist() %>% unique
  
  # do enrichment test here, and return list of enriched GO terms
  K = length(t.gns.neigh)
  
  #loop through every GO terms in the pop file
  if(K > 0){
    for (i in 1:nTerms) {
      GO.term.genes = all.GO.terms[i,which(all.GO.terms[i,] != '')]
      M = length(GO.term.genes)
      x.overlap.genes = intersect(GO.term.genes, t.gns.neigh)
      x = length(x.overlap.genes)
      info.gene.overrep[i,1] = rownames(all.GO.terms)[i]
      if(x > 0){
        info.gene.overrep[i,2] = phyper(x, M, N-M, K, lower.tail = FALSE) #He wrote : overlap.genes-1
        #insert FDR val
        info.gene.overrep[i,3] = p.adjust(info.gene.overrep[i,2], method = "bonferroni", n = nTerms)
      }
    }
    return(info.gene.overrep[which(info.gene.overrep[,2] < 0.05), 1])
  }else
    return("")
}
#----------- parallel excution
library(parallel)
numCores <- detectCores()
cl <- makeCluster(numCores)
clusterEvalQ(cl, {
  library(igraph)
  library(org.Hs.eg.db)
  library(stringi)
  library(AnnotationDbi)
  
})
clusterExport(cl=cl, varlist=c("all.GO.terms", "allTargetList", "getGOTerms", "getNeighbours", "V", "N", "ppi.net.all", "nTerms"), envir = .GlobalEnv)

result = parLapply(cl, seq_len(length(allTargetList)), fun = getGOTerms)
# result = parLapply(cl, seq_len(3), fun = getGOTerms)
stopCluster(cl)
saveRDS(result, file = "Data/result.MF.rds")

# hsGO <- godata('org.Hs.eg.db', ont="MF")
# findGOSim <- function(x, y, hsGO){
#   return(mgoSim(x, y, semData=hsGO, measure="Wang", combine="BMA"))
# }
# # pair-wise GO term similarity using GoSemSim
# library(foreach)
# library(doParallel)
# 
# cl <- makeCluster(detectCores())
# registerDoParallel(cl)
# 
# # clusterEvalQ(cl, {
# #   library(GOSemSim)
# # })
# # clusterExport(cl=cl, varlist=c("hsGO"), envir = .GlobalEnv)
# pairWiseGOSem.mat <- matrix(0, nrow = length(result), ncol = length(result))
# pairWiseGOSem.mat <- foreach(i = 1:length(result), .combine = 'rbind')  %:% 
#   foreach(j = 1:length(result), .combine = 'cbind', .packages = c('GOSemSim','foreach')) %dopar% {
#     pairWiseGOSem.mat[i,j] = findGOSim(result[[i]],result[[j]], hsGO)
#   }
# stopCluster(cl)
# 
# # pairWiseGOSem <- pblapply(result, function(x) lapply(result, findGOSim, x, hsGO))
# # pairWiseGOSem.mat <- do.call(rbind, lapply(pairWiseGOSem,function(x) do.call(cbind,x)))
# save(pairWiseGOSem.mat, file="pairWiseGOSem_MF.rda")
