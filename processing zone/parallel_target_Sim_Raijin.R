library(msa)
library(Biostrings)
library(seqinr)
library(foreach)
library(doParallel)

# ### Option list
args = commandArgs(trailingOnly=FALSE)

### Make Cluster nodes for parallelizing
dataDir = args[4]
dataFile = args[5]
startRow = strtoi(args[6])
offset = strtoi(args[7])
nClusters = strtoi(args[8])
outFile = args[9]

print(args)

cl <- makeCluster(nClusters)
registerDoParallel(cl)


## Read Fasta file
fasta <- readAAStringSet(paste0(dataDir, dataFile))
protSim <- function(seq1,seq2){
  return(Biostrings::pid(Biostrings::pairwiseAlignment(seq1,seq2)))
}

## 
endRow = startRow + offset -1
endRow = if (endRow > length(fasta)) length(fasta) else endRow
mat = matrix(0, nrow = (endRow - startRow + 1), ncol = length(fasta))


start_time <- Sys.time()
mat <- foreach(i = 1:nrow(mat), .combine='rbind') %:% 
  foreach(j = 1:ncol(mat), .combine='cbind', .packages='Biostrings') %dopar% {
    mat[i,j] = protSim(fasta[startRow + i - 1],fasta[j]);
  }
write.csv(mat, paste0(dataDir, outFile,"_",startRow, "_",endRow,".csv"), quote = F, row.names = F)
end_time <- Sys.time()
print(end_time - start_time)
stopCluster(cl)

