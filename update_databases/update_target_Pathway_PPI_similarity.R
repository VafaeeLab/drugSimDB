i2d.db = read.delim("../data/external/i2d.2_9.Public.HUMAN.tab", sep="\t", header = T)
drugRepo.db = read.delim("../data/new_net_info_V3_pval.csv", sep=",", header = T)
drugTarget.db = read.delim("../data/external/uniprot links.csv", sep=",", header = T)
drugTarget.db = drugTarget.db[which(drugTarget.db$Type == "SmallMoleculeDrug"),]

ppi.net.all = simplify(graph.data.frame(i2d.db[1:100,2:3], directed = F))

getNeighbours <- function(aNode){
  nb <- as.character(ego(g, order=1, nodes = aNode, mode = "all", mindist = 0)[[1]]$name)
  if(length(nb) == 1){# indicates no PPI parterns exists, except itself (since mindist=0)
    return(NA)
  }else{  # return the PPI partners except itself
    return(setdiff(nb,aNode))
  }
}