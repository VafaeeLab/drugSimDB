update_target_similarity <- function(outFile = "target_similarity"){
  
  require(msa)
  require(Biostrings)
  require(seqinr)
  require(pbapply)
  # require(foreach)
  # require(doParallel)
  # 
  # cl <- makeCluster(64)
  # registerDoParallel(cl)
  
  # Comute pairwise similarities between protein targets -----------------------
  
  # fasta <- readAAStringSet("./data/external/protein.fasta")
  # 
  # protSim <- function(seq1,seq2){
  #   return(Biostrings::pid(Biostrings::pairwiseAlignment(seq1,seq2)))
  # }
  # 
  # 
  # # prot.sim.matrix <- pbsapply(as(fasta, "list"), function(x) sapply(as(fasta, "list"), protSim, x))
  # 
  # # Multi-Thread version ----------------------------------------------------
  # 
  # library(parallel)
  # 
  # cl <- makeCluster(detectCores() - 1)
  # clusterEvalQ(cl, library(Biostrings))
  # clusterExport(cl, c("cl", "fasta", "protSim"), envir = .GlobalEnv)
  # 
  # prot.sim.matrix <- parallel::parSapply(cl, fasta, function(x) sapply(tmp, protSim, x))
  # 
  # stopCluster(cl)
  # 
  # ----------------------------------------------------
  
  
  prot.sim.matrix <- read.csv("./data/protein_similarity.csv")
  
  getIDs <- function(str){
    # return(gsub("drugbank_target\\|.","", strsplit(str," ")[[1]][1]))
    return(gsub("drugbank_target.","", str))
  }
  protIds <- sapply(as(colnames(prot.sim.matrix), "list"), getIDs)
  colnames(prot.sim.matrix) <- protIds
  rownames(prot.sim.matrix) <- protIds
  prot.sim.matrix = prot.sim.matrix/100
  
  # Read drug and their targets -------------------------------------------------------
  drug_target <- read.csv("./data/external/uniprot links.csv")
  drug_target <- subset(drug_target, Type == "SmallMoleculeDrug")[,c("DrugBank.ID", "UniProt.ID")]
  
  # compute target-based drug similarity -------------------------------------------------------
  
  # this function compute similarity between drug x and y based on their targets
  target.similarity <- function(x, y, drug_target, prot.sim.matrix){
    tx = as.character(drug_target[which(drug_target[,"DrugBank.ID"]==x),"UniProt.ID"])
    ty = as.character(drug_target[which(drug_target[,"DrugBank.ID"]==y),"UniProt.ID"])
    return(mean(prot.sim.matrix[tx,ty]))
  }
  
  tmp <- unique(subset(drug_target, select = "DrugBank.ID")[,1])
  ### following line (generating target similarity matrix) would be replaced by "newTargetSim.R 
  ### (+ newTargetSim.sh)" that runs on Raijin.
  ### Run this line only if the size of tmp is small enough
  # sim.matrix<- pbsapply(tmp, function(x) sapply(tmp, target.similarity, x, drug_target, prot.sim.matrix))
  
  ### Following line runs to combine batches of small chunks of the final matrix
  ### Ignore these lines (marked by START-END) if raijin isn't used 
  
  ### START
  library(gtools)
  datadir <- "C:\\Users\\Azad\\OneDrive - UNSW\\Vafaee Lab\\Projects\\Deep Brain\\Training 2\\targetSims\\"
  files <- mixedsort(list.files(path = datadir, pattern="*.csv", recursive=TRUE, full.names=TRUE), decreasing = F)
  
  sim.matrix <- NULL
  for(i in 1:length(files)){
    temp <- read.table(file = files[i], header = T, sep = ",")
    sim.matrix = rbind(sim.matrix,temp)
  }
  ### END
  
  
  #### run after all the batches finishes @ raijin 
  colnames(sim.matrix) <- tmp
  rownames(sim.matrix) <- tmp
  
  # Generate similarity matrix for all SmallMolecul drugs -------------------
  drugs.full <- read.csv("./data/external/drug links.csv", row.names = "DrugBank.ID")
  drugs.id.full <- rownames(subset(drugs.full, Drug.Type == "SmallMoleculeDrug"))
  
  full.sim.matrix <- matrix(data=NA,nrow=length(drugs.id.full),ncol=length(drugs.id.full))
  colnames(full.sim.matrix) <- drugs.id.full
  rownames(full.sim.matrix) <- drugs.id.full
  full.sim.matrix[rownames(sim.matrix), colnames(sim.matrix)] = as.matrix(sim.matrix)
  
  library(data.table)
  fwrite(full.sim.matrix, paste("./data/",outFile,".csv", sep =""), na = "NA", showProgress = T, verbose = T)
}






