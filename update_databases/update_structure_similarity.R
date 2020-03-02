update_structure_similarity <- function(outFile = "chem_similarity"){
  
  require(ChemmineR)
  require(fmcsR)
  require(pbapply)

  sdfset <- read.SDFset("./data/external/structures.sdf")
  valid <- validSDF(sdfset) # Identifies invalid SDFs in SDFset objects 
  sdfset <- sdfset[valid] # Removes invalid SDFs, if there are any
  drugbank_ids <- datablocktag(sdfset,1)
  
  
  apset <- sdf2ap(sdfset) # computes atom pair descriptors compounds
  valid <- which(sapply(as(apset, "list"), length)!=1)
  apset <- apset[valid] #view(apset[1:4])
  sdfset <- sdfset[valid] # Removes compounds failed to return APsfrom both apset and sdfset
  drugbank_ids <- datablocktag(sdfset,1)
  
  
  sim.matrix<- pbsapply(as(apset, "list"), function(x) sapply(as(apset, "list"), cmp.similarity, x))
  colnames(sim.matrix) <- drugbank_ids
  rownames(sim.matrix) <- drugbank_ids
 
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





