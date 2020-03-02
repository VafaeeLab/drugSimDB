main <- function(dataDir, FileTagList){
  require(parallel)
  
  for(i in 1:length(FileTagList)){
    load(file = paste0(dataDir,"pairWiseGoSem_",FileTagList[i]))  # loaded object name: pairWiseGOSem.mat
    
    numCores <- detectCores()
    cl <- makeCluster(numCores)
    
    # clusterEvalQ(cl, {
    #   library(igraph)
    # })
    
    clusterExport(cl=cl, varlist=c("pairWiseGOSem.mat"))
    
    result = parSapply()
    
  }
}