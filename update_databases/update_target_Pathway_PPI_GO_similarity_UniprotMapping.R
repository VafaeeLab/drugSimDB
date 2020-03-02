require(igraph)
require(GOSemSim)
require(pbapply)
require(clusterProfiler)
require(org.Hs.eg.db)
require(stringi)
# 1. Download GO terms from Gene Ontology website
goa.human = read.delim(file="data/external/goa_human_15thNov_2019.csv", sep="\t", header = F)
# As of July 2, 2019, only about 30% of all GO annotations were inferred computationally (Ref: http://www.geneontology.org/GO.consortiumlist.shtml)
# As these annotations are not checked by a human, the GO Consortium considers them to be marginally less reliable and they are commonly to higher
# level, less detailed terms (Ref: https://en.wikipedia.org/wiki/Gene_ontology#Terms_and_ontology)
# 2. Remove computationally inferred GO terms
goa.human = goa.human[-which(goa.human$V7 == 'IEA'),c(2,5,10)]
colnames(goa.human) = c("UniProt.ID", "GO.Term", "Description")
# split the db to get the GO term grouped for each protein: it becomes a list now
goa.human.list = split(as.character(goa.human[,2]), as.character(goa.human[,1]))
rm(goa.human)

# 1. Download GO terms from EnrichR library (https://amp.pharm.mssm.edu/Enrichr/#stats)
# BP
pop.filepath = "data/external/Enrichr/GO_Biological_Process_2018.txt"
num.col = max(count.fields(pop.filepath, sep = "\t", blank.lines.skip = TRUE), na.rm = TRUE)
GO_BP = read.table(pop.filepath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE, col.names = 1:num.col, blank.lines.skip = TRUE, quote = "")
# filter specific-general terms
nGS.terms.BP = rowSums(GO_BP != "", na.rm = T)-1
newBP = GO_BP[-which(nGS.terms.BP <= 15),]
newBP = newBP[-which(nGS.terms.BP >= 100),]
newBP = dplyr::select(newBP, -2) # throw away the the column with all NA (2nd-Column)

# MF
pop.filepath = "data/external/Enrichr/GO_Molecular_Function_2018.txt"
num.col = max(count.fields(pop.filepath, sep = "\t", blank.lines.skip = TRUE), na.rm = TRUE)
GO_MF = read.table(pop.filepath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE, col.names = 1:num.col, blank.lines.skip = TRUE, quote = "")
# filter specific-general terms
nGS.terms.MF = rowSums(GO_MF != "", na.rm = T)-1
newMF = GO_MF[-which(rowSums(GO_MF != "", na.rm = T)-1 <= 15),]
newMF = newMF[-which(rowSums(newMF != "", na.rm = T)-1 >= 100),]
newMF = dplyr::select(newMF, -2) # throw away the the column with all NA (2nd-Column)

# CC
pop.filepath = "data/external/Enrichr/GO_Cellular_Component_2018.txt"
num.col = max(count.fields(pop.filepath, sep = "\t", blank.lines.skip = TRUE), na.rm = TRUE)
GO_CC = read.table(pop.filepath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, fill = TRUE, col.names = 1:num.col, blank.lines.skip = TRUE, quote = "")
newCC = GO_CC[-which(rowSums(GO_CC != "", na.rm = T)-1 <= 15),]
newCC = newCC[-which(rowSums(newCC != "", na.rm = T)-1 >= 100),]
newCC = dplyr::select(newCC, -2) # throw away the the column with all NA (2nd-Column)

# 2. Combine BP, MF, CC
all.GO.terms = bind_rows(newBP,newMF,newCC)
# Parse the GO accession numbers for each Terms
row.names(all.GO.terms) = gsub(")","",lapply(lapply(stri_split_fixed(stri_reverse(all.GO.terms[,1]),"(",n = 2), FUN = stri_reverse), "[[",1))
nTerms = nrow(all.GO.terms)
pop.genes = unique(c(as.matrix(all.GO.terms[,-1])))
N = length(pop.genes)

# 
# Download Human PPI network from i2d database
i2d.db = read.delim("data/external/i2d.2_9.Public.HUMAN.tab", sep="\t", header = T)
# Load the PPI as an iGraph object and remove cycles and loops
ppi.net.all = igraph::simplify(graph.data.frame(i2d.db[,c("SwissProt1","SwissProt2")], directed = F))
rm(i2d.db)

# The drug-pairs that matters
drugRepo.db = read.delim("data/new_net_info_V3_pval.csv", sep=",", header = T)
# The target proteins of each Drugs
drugTarget.db = read.delim("data/external/uniprot links.csv", sep=",", header = T)
# only focus on the Small Molecule drugs
drugTarget.db = drugTarget.db[which(drugTarget.db$Type == "SmallMoleculeDrug"),c(1,4)]
# split the db to get the target grouped for each drugs: it becomes a list now
drugTarget.list = split(as.character(drugTarget.db[,2]), as.character(drugTarget.db[,1]))
rm(drugTarget.db)

# Get PPI neighbours of a Protein
getNeighbours <- function(aNode){
  #if the node doesn't exist in the net, return null
  if(sum(which(V(ppi.net.all)$name == aNode)) == 0)
    return(NA)
  nb <- as.character(ego(ppi.net.all, order=1, nodes = aNode, mode = "all", mindist = 0)[[1]]$name)
  if(length(nb) == 1){# indicates no PPI parterns exists, except itself (since mindist=0)
    return(NA)
  }else{  # return the PPI partners except itself
    # return(setdiff(nb,aNode))
    return (nb)
  }
}
getGOTerms <- function(aDrug){
  # get all the targets associated with the aDrug
  t.prots = aDrug[]
  # print(t.prots)
  # find a Set of PPI neighbours of all the target proteins
  t.prots.neigh = unique(unlist(lapply(t.prots, FUN = function(x){return(getNeighbours(x))})))
  # get list of GO terms associated with these proteins
  t.GOs = lapply(t.prots.neigh, FUN = function(x){
    return(goa.human.list[x])
  }) %>% base::unlist() %>% unique()
  
  # print(GOterms)
  return(t.GOs)
}

result.CC = pblapply(drugTarget.list[["DB00116"]], FUN = getGOTerms)
# result = lapply(drugTarget.list["DB00194"], FUN = getGOTerms)
# hsGO_BP <- godata('org.Hs.eg.db', ont="BP")
hsGO_CC <- godata('org.Hs.eg.db', ont="CC")
# hsGO_MF <- godata('org.Hs.eg.db', ont="MF")
findGOSim <- function(x, y, hsGO){
  return(mgoSim(x[], y[], semData=hsGO, measure="Wang", combine="BMA"))
}
# # pair-wise GO term similarity using GoSemSim
pairWiseGOSem_CC <- lapply(result.CC, function(x) pblapply(result.CC, findGOSim, x, hsGO_CC))
# # ## testing
# nodes_of_interest = c("P63104", "Q9H2F3", "A0AV96")
# 
# 
