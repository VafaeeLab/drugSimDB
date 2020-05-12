library("parallel")
library("foreach")
library("doParallel")

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
JI.similarity <- function(x,y){
  jaccard = NA
  if(length(union(x,y))!=0){jaccard = length(intersect(x,y))/length(union(x,y))}
  return(jaccard)
}

## 
endRow = startRow + offset -1
endRow = if (endRow > length(result)) length(result) else endRow
mat = matrix(0, nrow = (endRow - startRow + 1), ncol = length(result))


start_time <- Sys.time()

mat <- foreach(i = 1:nrow(mat), .combine = 'rbind')  %:% 
  foreach(j = 1:ncol(mat), .combine = 'cbind', .packages = c('foreach')) %dopar% {
        mat[i,j] = JI.similarity(result[[startRow + i - 1]],result[[j]])
  }

# Save the partial matrix
write.csv(mat, paste0(dataDir, outFile,"_",ontTag,"_",startRow, "_",endRow,".csv"), quote = F, row.names = F)
end_time <- Sys.time()
print(end_time - start_time)
stopCluster(cl)

