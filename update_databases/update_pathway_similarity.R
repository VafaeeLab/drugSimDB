update_pathway_similarity <- function(outFile = "path_similarity"){
  
  # KEGG pathways associated with each drug: Uncomment if update is needed -------------------------------------------
  # source("./update_KEGGpathways_for_drugs.R")
  # update_KEGGpathways_for_drugs()
  
  # Compute pathway similarities ------------------------------------------------------------------------
  
  require(pbapply)
  drug_path <- read.csv("./data/durg_KEGGpathways.csv")
  
  path.similarity <- function(x,y, drug_path){
    px = drug_path[which(drug_path[,"Drug.ID"]==x),"Kegg.ID"]
    py = drug_path[which(drug_path[,"Drug.ID"]==y),"Kegg.ID"]
    jaccard = NA
    if(length(union(px,py))!=0){jaccard = length(intersect(px,py))/length(union(px,py))}
    return(jaccard)
  }
  
  tmp <- unique(subset(drug_path, select = "Drug.ID")[,1])
  sim.matrix<- pbsapply(tmp, function(x) sapply(tmp, path.similarity, x, drug_path))
  
  colnames(sim.matrix) <- tmp
  rownames(sim.matrix) <- tmp
  
  # Generate similarity matrix for all SmallMolecul drugs -------------------
  drugs.full <- read.csv("./data/external/drug links.csv", row.names = "DrugBank.ID")
  drugs.id.full <- rownames(subset(drugs.full, Drug.Type == "SmallMoleculeDrug"))
  
  full.sim.matrix <- matrix(data=NA,nrow=length(drugs.id.full),ncol=length(drugs.id.full))
  colnames(full.sim.matrix) <- drugs.id.full
  rownames(full.sim.matrix) <- drugs.id.full
  full.sim.matrix[rownames(sim.matrix), colnames(sim.matrix)] = sim.matrix
  
  # Save the similarity matrix ----------------------------------------------
  
  write.csv(full.sim.matrix, paste("./data/",outFile,".csv", sep =""))
  
}





