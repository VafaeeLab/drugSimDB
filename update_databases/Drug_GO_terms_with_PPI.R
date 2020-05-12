require(igraph)
require(GOSemSim)
require(pbapply)
require(clusterProfiler)
require(org.Hs.eg.db)
require(stringi)
require(AnnotationDbi)
library(foreach)
library(doParallel)
library(parallel)

setwd("/srv/scratch/z3526914/DrugRepo/")
args = commandArgs(trailingOnly=FALSE)
ontTag <- args[6]
pop.filepath = args[7]
cutoff = as.numeric(args[8])


# 1. Download GO terms from EnrichR library (https://amp.pharm.mssm.edu/Enrichr/#stats)

num.col = max(count.fields(pop.filepath, sep = "\t", blank.lines.skip = TRUE), na.rm = TRUE)
GO = read.table(pop.filepath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE, col.names = 1:num.col, blank.lines.skip = TRUE, quote = "")
# filter specific-general terms
nGS.terms = rowSums(GO != "", na.rm = T)-1
new = GO[-which(rowSums(GO != "", na.rm = T)-1 <= 15),]
new = new[-which(rowSums(new != "", na.rm = T)-1 >= 100),]
new = dplyr::select(new, -2) # throw away the the column with all NA (2nd-Column)

# 2. 
all.GO.terms = new
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
    return(info.gene.overrep[which(info.gene.overrep[,2] < cutoff), 1])
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
clusterExport(cl=cl, varlist=c("all.GO.terms", "allTargetList", "getGOTerms", "getNeighbours", "V", "N", "ppi.net.all", "nTerms", "cutoff"), envir = .GlobalEnv)
result = parLapply(cl, seq_len(length(allTargetList)), fun = getGOTerms)
stopCluster(cl)
saveRDS(result, file = paste0("Data/result.", ontTag ,".JI.with.PPI.rds"))


# cl <- makeCluster(numCores)
# JI.similarity <- function(x,y){
#   jaccard = NA
#   if(length(union(x,y))!=0){jaccard = length(intersect(x,y))/length(union(x,y))}
#   return(jaccard)
# }
# mat = matrix(0, nrow = length(result), ncol = length(result))
# start_time <- Sys.time()
# 
# mat <- foreach(i = 1:nrow(mat), .combine = 'rbind')  %:% 
#   foreach(j = 1:ncol(mat), .combine = 'cbind', .packages = c('foreach')) %dopar% {
#     mat[i,j] = JI.similarity(result[[i]],result[[j]])
#   }
# 
# # Save the partial matrix
# write.csv(mat, paste0("/srv/scratch/z3526914/DrugRepo/Data/GO_sim_JI_with_PPI_",ontTag,".csv"), quote = F, row.names = F)
# end_time <- Sys.time()
# print(end_time - start_time)
# stopCluster(cl)

