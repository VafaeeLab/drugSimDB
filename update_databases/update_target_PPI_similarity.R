update_target_PPI_similarity <- function(outFile = "target_PPI_similarity"){
  
  require(pbapply)
  require(doParallel)
  require(doSNOW)
  # Read drug and their targets -------------------------------------------------------
  uni.db <- read.csv("./data/external/uniprot links.csv")
  uni.db <- subset(uni.db, Type == "SmallMoleculeDrug")[,c("DrugBank.ID", "UniProt.ID")]
  drugs <- unique(uni.db$DrugBank.ID)
  
  # read Human PPI from i2d database
  i2d.db = read.delim("./data/external/i2d.2_9.Public.HUMAN.tab", header = TRUE, sep = "\t")
  
  # get PPI info
  getPPI <- function(m,n){
    count = 0
    for(i in 1:length(m)){
      print(m[i])
      for(j in 1:length(n)){
        print(n[j])
        if((as.character(i2d.db$SwissProt1) == m[i] & as.character(i2d.db$SwissProt2) == n[j]) | 
           (as.character(i2d.db$SwissProt1) == n[j] & as.character(i2d.db$SwissProt2) == m[i]))
          count = count + 1
      }
    }
  }
  # find similarity
  mat <- matrix(0, nrow = length(drugs), ncol = length(drugs))
  
  target.ppi.similarity <- function(x,y){
    set_x = as.character(uni.db[which(as.character(uni.db$DrugBank.ID) == x),"UniProt.ID"])
    set_y = as.character(uni.db[which(as.character(uni.db$DrugBank.ID) == y),"UniProt.ID"])
    # print(set_x)
    # print(set_y)
    # set_x = setdiff(set_x,set_y)
    # set_y = setdiff(set_y,set_x)
    x_PPIs = i2d.db[which(as.character(i2d.db$SwissProt1) %in% set_x),c(2,3)]
    y_PPIs = i2d.db[which(as.character(i2d.db$SwissProt2) %in% set_y),c(2,3)]
    
    x_val_1 = apply(x_PPIs, 1, FUN = function(aRow){
      return(paste0(aRow[1], ",", aRow[2]))
    })
    y_val_1 = apply(y_PPIs, 1, FUN = function(aRow){
      return(paste0(aRow[1], ",", aRow[2]))
    })
    x_val_2 = apply(x_PPIs, 1, FUN = function(aRow){
      return(paste0(aRow[2], ",", aRow[1]))
    })
    y_val_2 = apply(y_PPIs, 1, FUN = function(aRow){
      return(paste0(aRow[2], ",", aRow[1]))
    })
    
    nPPI <- length(intersect(x_val_1,y_val_1)) + length(intersect(x_val_2,y_val_2))
    return(nPPI/(length(set_x)*length(set_y)))
  }
  
  # try 
  mat <- foreach(i = 1:nrow(mat), .combine='rbind') %:% 
    foreach(j = 1:ncol(mat), .combine='cbind') %dopar% {
      if(i != j)
        {mat[i,j] = target.ppi.similarity(drugs[i],drugs[j])}
    }
  save(mat,file = "PPI_mat.rda")
}