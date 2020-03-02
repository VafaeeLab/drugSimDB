library(foreach)
library(doParallel)
library(parallel)
library(GOSemSim)

# ### Option list
args = commandArgs(trailingOnly=FALSE)
dataDir = args[6]
dataFile = args[7]
startRow = strtoi(args[8])
offset = strtoi(args[9])
ontTag = args[10]
nClusters = strtoi(args[11])
outFile = args[12]

print(args)

cl <- makeCluster(nClusters)
registerDoParallel(cl)

# read data
result <- readRDS(paste0(dataDir, dataFile))
hsGO <- godata('org.Hs.eg.db', ont=ontTag)
findGOSim <- function(x, y, hsGO){
  return(mgoSim(x, y, semData=hsGO, measure="Wang", combine="BMA"))
}

## 
endRow = startRow + offset -1
endRow = if (endRow > length(result)) length(result) else endRow
mat = matrix(0, nrow = (endRow - startRow + 1), ncol = length(result))


start_time <- Sys.time()

mat <- foreach(i = 1:nrow(mat), .combine = 'rbind')  %:% 
  foreach(j = 1:ncol(mat), .combine = 'cbind', .packages = c('GOSemSim','foreach')) %dopar% {
    mat[i,j] = findGOSim(result[[startRow + i - 1]],result[[j]], hsGO)
  }

# Save the partial matrix
write.csv(mat, paste0(dataDir, outFile,"_",ontTag,"_",startRow, "_",endRow,".csv"), quote = F, row.names = F)
end_time <- Sys.time()
print(end_time - start_time)
stopCluster(cl)

